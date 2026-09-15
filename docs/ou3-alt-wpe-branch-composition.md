# ALT WPE branch composition

The startup and Live products carry the exact WPE shadow and two machine WPE
histories independently. Log initialization, valid-period production and the
usable-period latch need not agree. The machine state contains both moment
histories, both log histories and both latches. Its source-owned predecessor
and successor are retained across every represented event.

## Sample-entry frequency

For each machine history, `MachineRead` represents the literal wrapper order:

1. Read that history's stored log. If finite, evaluate the frequency getter at
   its negation, even when the usable latch is false; otherwise return NaN.
2. Select the getter result only if the same history's latch is true and the
   returned value is finite and positive. Otherwise select the compiled prior.
3. Feed this selected value into the existing machine frontend and tuner.

The exact shadow selects its own frequency, which must equal the raw getter
input to the same lower event; equality after clamping is insufficient. Its
selected value is never used to decide the machine branch. Their exact difference remains in the supply
ledger. A nonfinite getter result is an explicit return class; its selection
falls back to the prior even after the latch has fired. The getter and pre-clamp
input accept the full finite binary32 lattice,
including positive subnormals; the statistics floor precedes normal arithmetic.
Nonfinite log production and arithmetic outside the named finite graph are not
thereby qualified.

The downstream consumers retain shipping order and distinct operands:

| Consumer | Frequency source |
| --- | --- |
| Racc | Current sample-entry machine getter/prior, before the two clamps |
| Adaptive band | Previous machine statistics frequency, with the external fallback and band bounds |
| Current statistics | Current machine getter/prior, through statistics bounds |
| Tau and subsequent tuner candidate | Current stored statistics frequency through the outer tuning bounds |

`finite_machine_frontend_sigma_source.bind_step` checks these operands against
the carried frontend and the same updated band/statistics state. The Racc join
reads the raw selected frequency from the same frequency result. Changing one
mode's takeover choice therefore changes all of that mode's downstream inputs;
it cannot retain a frontend computed for the other choice.

## Current-sample period production

The lower exact physical/filter event executes once. Its already-produced
machine Mahony vertical sample drives both machine WPE moment histories. Each
history determines its own raw-period branch and its own log successor. On a
valid period, the log input must be the raw period from those same machine
moments and square-root result. The post-update usability getter is attached
to the updated log in the same history.

The lower tuner/log projection and the full moment/log/latch product are joined
by equality of their machine log predecessors and successors. Mode-specific
log witnesses in the lower projection do not independently authorize a branch:
the full moment step must produce that branch and the same successor. A held
machine branch consumes no log-update witness even when the exact shadow
produces a period. Conversely, a machine branch can produce a period while the
exact shadow holds. The two machine histories can also initialize on different
samples.

This relation has no evaluation cycle: the machine tuner reads only sample-entry
WPE state, and the current WPE update follows the tuner. The wrapper's evaluation
order checks one conjunction of relations, without replaying the physical or
filter event or replacing its source packet.

## Prefix composition

Construction installs the literal WPE reset; the startup wrapper rejects a
replacement full-WPE root. An IMU transition requires the same persistent
predecessor, the same machine dt and the same Mahony vertical sample, then binds
the full successor into the next state. goLive, MAG and HOLD preserve that
machine state by identity. Induction therefore preserves these relations for
every represented finite prefix. The existing complete-word constructor requires
600 additional full-WPE IMU transitions after the startup-produced Live entry.

Exact/machine comparison agreement is thus not a premise of this product. The
old branch-robustness obligation is replaced by the composed independent branch
relation, not by an assertion that rounding preserves every comparison.

## Qualification boundary

This closes branch composition on the represented arithmetic graph. It does not
prove that every admitted source executes wholly inside that graph, that startup
reaches Live by its conditional deadline, or that the target ESP32/Eigen/libm
implements the named arithmetic relations. Uniform arithmetic supplies, all
exceptional outcomes and the complete source-uniform 600-step word remain
separate master prerequisites. No rho or ALT PASS follows from this relation.
