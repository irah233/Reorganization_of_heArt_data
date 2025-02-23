from pyadjoint.drivers import compute_gradient, compute_hessian
from pyadjoint.enlisting import Enlist
from pyadjoint.tape import get_working_tape, stop_annotating, no_annotations
import numpy as np
from dolfin import MPI 
from dolfin import HDF5File


class ReducedFunctional(object):
    """Class representing the reduced functional.

    A reduced functional maps a control value to the provided functional.
    It may also be used to compute the derivative of the functional with
    respect to the control.

    Args:
        functional (:obj:`OverloadedType`): An instance of an OverloadedType,
            usually :class:`AdjFloat`. This should be the return value of the
            functional you want to reduce.
        controls (list[Control]): A list of Control instances, which you want
            to map to the functional. It is also possible to supply a single Control
            instance instead of a list.

    """

    def __init__(self, functional, controls, scale=1.0, tape=None,
                 eval_cb_pre=lambda *args: None,
                 eval_cb_post=lambda *args: None,
                 derivative_cb_pre=lambda *args: None,
                 derivative_cb_post=lambda *args: None,
                 hessian_cb_pre=lambda *args: None,
                 hessian_cb_post=lambda *args: None,
                 debugfile=None):

        self.functional = functional
        self.tape = get_working_tape() if tape is None else tape
        self.controls = Enlist(controls)
        self.scale = scale
        self.eval_cb_pre = eval_cb_pre
        self.eval_cb_post = eval_cb_post
        self.derivative_cb_pre = derivative_cb_pre
        self.derivative_cb_post = derivative_cb_post
        self.hessian_cb_pre = hessian_cb_pre
        self.hessian_cb_post = hessian_cb_post
        self.debugfile = debugfile

    def derivative(self, options={}):
        """Returns the derivative of the functional w.r.t. the control.

        Using the adjoint method, the derivative of the functional with
        respect to the control, around the last supplied value of the control,
        is computed and returned.

        Args:
            options (dict): A dictionary of options. To find a list of available options
                have a look at the specific control type.

        Returns:
            OverloadedType: The derivative with respect to the control.
                Should be an instance of the same type as the control.

        """
        # Call callback
        values = [c.data() for c in self.controls]
        self.derivative_cb_pre(self.controls.delist(values))

        derivatives = compute_gradient(self.functional,
                                       self.controls,
                                       options=options,
                                       tape=self.tape,
                                       adj_value=self.scale)

        # Call callback
        self.derivative_cb_post(self.functional.block_variable.checkpoint,
                                self.controls.delist(derivatives),
                                self.controls.delist(values))

        return self.controls.delist(derivatives)

    @no_annotations
    def hessian(self, m_dot, options={}):
        """Returns the action of the Hessian of the functional w.r.t. the control on a vector m_dot.

        Using the second-order adjoint method, the action of the Hessian of the
        functional with respect to the control, around the last supplied value
        of the control, is computed and returned.

        Args:
            m_dot ([OverloadedType]): The direction in which to compute the
                action of the Hessian.
            options (dict): A dictionary of options. To find a list of
                available options have a look at the specific control type.

        Returns:
            OverloadedType: The action of the Hessian in the direction m_dot.
                Should be an instance of the same type as the control.
        """
        # Call callback
        values = [c.data() for c in self.controls]
        self.hessian_cb_pre(self.controls.delist(values))

        r = compute_hessian(self.functional, self.controls, m_dot, options=options, tape=self.tape)

        # Call callback
        self.hessian_cb_post(self.functional.block_variable.checkpoint,
                             self.controls.delist(r),
                             self.controls.delist(values))

        return self.controls.delist(r)

    @no_annotations
    def __call__(self, values):
        """Computes the reduced functional with supplied control value.

        Args:
            values ([OverloadedType]): If you have multiple controls this should be a list of
                new values for each control in the order you listed the controls to the constructor.
                If you have a single control it can either be a list or a single object.
                Each new value should have the same type as the corresponding control.

        Returns:
            :obj:`OverloadedType`: The computed value. Typically of instance
                of :class:`AdjFloat`.

        """
        values = Enlist(values)
        if len(values) != len(self.controls):
            raise ValueError("values should be a list of same length as controls.")

        # Call callback.
        self.eval_cb_pre(self.controls.delist(values))

        comm = values[0].function_space().mesh().mpi_comm()
        rank = comm.Get_rank()
        fdata = HDF5File(comm, self.debugfile, 'w')

        for i, value in enumerate(values):
            self.controls[i].update(value)
            value_ = self.GetValue(value)
            fdata.write(value, "Tmax"+str(i)) 

        fdata.close()

        self.tape.reset_blocks()
        blocks = self.tape.get_blocks()
        with self.marked_controls():
            with stop_annotating():
                for i in range(len(blocks)):
                    blocks[i].recompute()

        func_value = self.scale * self.functional.block_variable.checkpoint

        # Call callback
        self.eval_cb_post(func_value, self.controls.delist(values))

        return func_value

    def optimize_tape(self):
        self.tape.optimize(
            controls=self.controls,
            functionals=[self.functional]
        )

    def marked_controls(self):
        return marked_controls(self)

    def GetValue(self, var):

        comm = var.function_space().mesh().mpi_comm()

        try:
            val_local = var.vector().get_local()[0]
        except IndexError:
            val_local = 0.0
        
        
        value = MPI.sum(comm, val_local)
        
        return value




class marked_controls(object):
    def __init__(self, rf):
        self.rf = rf

    def __enter__(self):
        for control in self.rf.controls:
            control.mark_as_control()

    def __exit__(self, *args):
        for control in self.rf.controls:
            control.unmark_as_control()
