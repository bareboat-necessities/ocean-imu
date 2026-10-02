"""Reachable-graph contract for source-uniform OU-III shaped storage.

This is a fail-closed theorem schema. Generated outputs may not be supplied as
independent boxes; a certificate must propagate one shared history cell.
"""
from __future__ import annotations
from dataclasses import dataclass

GENERATED=("mahony","frequency","variance","tau","sigma_aw","R_S","T_S",
           "covariance","gain","scheduler_phase","mag_reference","ba_gate")

@dataclass(frozen=True)
class ReachableGraphCertificate:
    horizon_s: float
    shared_history_input: bool
    generated_outputs_propagated_causally: bool
    persistent_cross_word_state: bool
    homogeneous_zero_set_excluded: bool
    compact_reachable_graph: bool
    linked_supply_enclosed: bool
    every_prefix_enclosed: bool
    entry_comparison_to_local_V: bool

    def validates_source_uniform_existence(self):
        return (self.shared_history_input and self.generated_outputs_propagated_causally
                and self.persistent_cross_word_state and self.homogeneous_zero_set_excluded
                and self.compact_reachable_graph)

    def validates_constructive_entry(self):
        return (self.validates_source_uniform_existence() and self.linked_supply_enclosed
                and self.every_prefix_enclosed and self.entry_comparison_to_local_V)

def reject_cartesian_generated_boxes(spec: dict):
    bad=[k for k in GENERATED if k in spec and isinstance(spec[k],(tuple,list))
         and len(spec[k])==2]
    if bad:
        raise ValueError("independent generated-output boxes forbidden: "+",".join(bad))
    return True
