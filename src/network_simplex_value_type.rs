use ebi_arithmetic::{Fraction, fraction::fraction_f64::FractionF64};
use malachite::{Integer, Natural};

pub trait MulWithFloat {
    fn mul_with_float(self, rhs: &f64) -> Self;
}

impl MulWithFloat for f64 {
    fn mul_with_float(self, rhs: &f64) -> Self {
        self * rhs
    }
}

impl MulWithFloat for i64 {
    fn mul_with_float(self, _rhs: &f64) -> Self {
        // this should never occur. it is necessary to make network simplex work on both integers and floats
        panic!("Cannot multiply values of different types");
    }
}

impl MulWithFloat for i128 {
    fn mul_with_float(self, _rhs: &f64) -> Self {
        // this should never occur. it is necessary to make network simplex work on both integers and floats
        panic!("Cannot multiply values of different types");
    }
}

impl MulWithFloat for Integer {
    fn mul_with_float(self, _rhs: &f64) -> Self {
        // this should never occur. it is necessary to make network simplex work on both integers and floats
        panic!("Cannot multiply values of different types");
    }
}

impl MulWithFloat for Natural {
    fn mul_with_float(self, _rhs: &f64) -> Self {
        // this should never occur. it is necessary to make network simplex work on both integers and floats
        panic!("Cannot multiply values of different types");
    }
}

impl MulWithFloat for FractionF64 {
    fn mul_with_float(self, rhs: &f64) -> Self {
        self * *rhs
    }
}

pub trait ToBigInt {
    fn to_big_int(&self) -> Integer;
}

impl ToBigInt for f64 {
    // this should never occur. it is necessary to make network simplex work on both integers and floats
    fn to_big_int(&self) -> Integer {
        panic!("Cannot multiply values of different types");
    }
}

impl ToBigInt for i64 {
    fn to_big_int(&self) -> Integer {
        Integer::from(*self)
    }
}

impl ToBigInt for i128 {
    fn to_big_int(&self) -> Integer {
        Integer::from(*self)
    }
}

impl ToBigInt for Integer {
    fn to_big_int(&self) -> Integer {
        self.clone()
    }
}

impl ToBigInt for Fraction {
    fn to_big_int(&self) -> Integer {
        panic!("Cannot multiply values of different types");
    }
}
