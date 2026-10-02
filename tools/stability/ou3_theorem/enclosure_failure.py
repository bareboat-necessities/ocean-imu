"""Typed fail-closed exception for causal proof enclosures."""
class EnclosureFailure(ArithmeticError):
 def __init__(self,stage,reason,time=None,input_sensitivity=None):
  self.stage=str(stage);self.reason=str(reason);self.time=time
  self.input_sensitivity=input_sensitivity or {}
  msg=f"{self.stage} enclosure unresolved"
  if time is not None:msg+=f" at t={float(time):.6f}"
  msg+=f": {self.reason}"
  super().__init__(msg)
