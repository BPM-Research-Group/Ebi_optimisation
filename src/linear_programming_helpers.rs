use ebi_arithmetic::Zero;
use sprs::{CsVecBase, CsVecView};
use std::ops::Deref;

pub(crate) fn resized_view<IStorage, DStorage>(
    vec: &CsVecBase<IStorage, DStorage, AbnormalFraction>,
    len: usize,
) -> CsVecView<AbnormalFraction>
where
    IStorage: Deref<Target = [usize]>,
    DStorage: Deref<Target = [AbnormalFraction]>,
{
    let mut indices = vec.indices();
    let mut data = vec.data();
    while let Some(&i) = indices.last() {
        if i < len {
            // TODO: binary search
            break;
        }

        indices = &indices[..(indices.len() - 1)];
        data = &data[..(data.len() - 1)];
    }

    // Safety: new indices and data are the same size,indices are still sorted and all indices
    // are less than the new length. Thus, all CsVecView invariants are satisfied.
    unsafe { CsVecView::new_uncheked(len, indices, data) }
    // unsafe { CsVecView::new_view_raw(len, data.len(), indices.as_ptr(), data.as_ptr()) }
}

pub(crate) fn to_dense<IStorage, DStorage>(
    vec: &CsVecBase<IStorage, DStorage, AbnormalFraction>,
) -> Vec<AbnormalFraction>
where
    IStorage: Deref<Target = [usize]>,
    DStorage: Deref<Target = [AbnormalFraction]>,
{
    let mut dense = vec![f0_ab!(); vec.dim()];
    vec.scatter(&mut dense);
    dense
}

#[cfg(test)]
use sprs::{CsMat, CsVec};

use crate::{abnormal_fraction::AbnormalFraction, f0_ab};

#[cfg(test)]
pub(crate) fn to_sparse(slice: &[AbnormalFraction]) -> CsVec<AbnormalFraction> {
    let mut res = CsVec::empty(slice.len());
    for (i, val) in slice.iter().enumerate() {
        use ebi_arithmetic::Zero;

        if !val.is_zero() {
            res.append(i, val.clone());
        }
    }
    res
}

#[cfg(test)]
pub(crate) fn assert_matrix_eq(mat: &CsMat<AbnormalFraction>, reference: &[Vec<AbnormalFraction>]) {
    let mat = mat.to_csr();
    assert_eq!(mat.rows(), reference.len());
    for (r, row) in mat.outer_iterator().enumerate() {
        assert_eq!(to_dense(&row), reference[r], "matrices differ in row {}", r);
    }
}
