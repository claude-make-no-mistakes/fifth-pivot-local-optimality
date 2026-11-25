import Init

namespace Pivot

/-! The relaxed-family necessary condition from the fourth-pivot cap. -/

theorem necessary_product (bench growth fourth ratio cap : Rat)
    (gain : bench < growth)
    (stage : growth = (1 + ratio) * fourth)
    (nonnegative : 0 ≤ 1 + ratio)
    (fourth_bound : fourth ≤ cap) :
    bench < (1 + ratio) * cap := by
  have comparison := Rat.mul_le_mul_of_nonneg_left fourth_bound nonnegative
  rw [← stage] at comparison
  exact Rat.not_le.mp fun le => Rat.not_le.mpr gain (Rat.le_trans comparison le)

theorem local_strict (bench growth comparator : Rat)
    (increasing : growth < comparator)
    (bound : comparator ≤ bench) :
    growth < bench :=
  Rat.not_le.mp fun le => Rat.not_le.mpr increasing (Rat.le_trans bound le)

end Pivot

#print axioms Pivot.necessary_product
#print axioms Pivot.local_strict
