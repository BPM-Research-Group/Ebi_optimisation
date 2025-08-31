use std::{
    cmp::Ordering,
    fmt::Display,
    iter::Sum,
    ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign},
};

use anyhow::anyhow;
use ebi_arithmetic::{Fraction, MaybeExact, One, Round, Signed, Zero};

#[derive(Eq, PartialEq, Clone, Debug)]
pub enum AbnormalFraction {
    Normal(Fraction),
    Infinite,
    NegInfinite,
    NaN,
}

impl AbnormalFraction {
    pub fn infinity() -> Self {
        Self::Infinite
    }

    pub fn neg_infinity() -> Self {
        Self::NegInfinite
    }

    /// Returns `true` if this number is neither infinite nor NaN.
    pub fn is_finite(&self) -> bool {
        match self {
            AbnormalFraction::Normal(_) => true,
            AbnormalFraction::Infinite | AbnormalFraction::NegInfinite => false,
            AbnormalFraction::NaN => false,
        }
    }

    /// Returns true if this value is positive infinity or negative infinity, and false otherwise.
    pub fn is_infinite(&self) -> bool {
        match self {
            AbnormalFraction::Normal(_) => false,
            AbnormalFraction::Infinite | AbnormalFraction::NegInfinite => true,
            AbnormalFraction::NaN => false,
        }
    }

    pub fn is_neg_infinite(&self) -> bool {
        match self {
            AbnormalFraction::Normal(_) | AbnormalFraction::Infinite | AbnormalFraction::NaN => {
                false
            }
            AbnormalFraction::NegInfinite => true,
        }
    }

    pub(crate) fn both_normal(&self, rhs: &Self) -> bool {
        match (self, rhs) {
            (AbnormalFraction::Normal(_), AbnormalFraction::Normal(_)) => true,
            _ => false,
        }
    }
}

impl Display for AbnormalFraction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AbnormalFraction::Normal(fr) => fr.fmt(f),
            AbnormalFraction::Infinite => write!(f, "∞"),
            AbnormalFraction::NegInfinite => write!(f, "-∞"),
            AbnormalFraction::NaN => write!(f, "NaN"),
        }
    }
}

impl Signed for AbnormalFraction {
    fn abs(self) -> Self {
        match self {
            AbnormalFraction::Normal(f) => AbnormalFraction::Normal(f.abs()),
            AbnormalFraction::Infinite => Self::NegInfinite,
            AbnormalFraction::NegInfinite => Self::Infinite,
            AbnormalFraction::NaN => Self::NaN,
        }
    }

    fn is_positive(&self) -> bool {
        match self {
            AbnormalFraction::Normal(f) => f.is_positive(),
            AbnormalFraction::Infinite => true,
            AbnormalFraction::NegInfinite => false,
            AbnormalFraction::NaN => false,
        }
    }

    fn is_negative(&self) -> bool {
        match self {
            AbnormalFraction::Normal(f) => f.is_negative(),
            AbnormalFraction::Infinite => false,
            AbnormalFraction::NegInfinite => true,
            AbnormalFraction::NaN => false,
        }
    }

    fn is_not_negative(&self) -> bool {
        match self {
            AbnormalFraction::Normal(f) => f.is_not_negative(),
            AbnormalFraction::Infinite => true,
            AbnormalFraction::NegInfinite => false,
            AbnormalFraction::NaN => true,
        }
    }

    fn is_not_positive(&self) -> bool {
        match self {
            AbnormalFraction::Normal(f) => f.is_not_positive(),
            AbnormalFraction::Infinite => false,
            AbnormalFraction::NegInfinite => true,
            AbnormalFraction::NaN => true,
        }
    }
}

impl Default for AbnormalFraction {
    fn default() -> Self {
        Self::Normal(Fraction::zero())
    }
}

impl Zero for AbnormalFraction {
    fn zero() -> Self {
        Self::Normal(Fraction::zero())
    }

    fn is_zero(&self) -> bool {
        match self {
            AbnormalFraction::Normal(f) => f.is_zero(),
            AbnormalFraction::Infinite | AbnormalFraction::NegInfinite | AbnormalFraction::NaN => {
                false
            }
        }
    }
}

impl num_traits::identities::Zero for AbnormalFraction {
    fn zero() -> Self {
        Self::Normal(Fraction::zero())
    }

    fn is_zero(&self) -> bool {
        match self {
            AbnormalFraction::Normal(f) => f.is_zero(),
            AbnormalFraction::Infinite | AbnormalFraction::NegInfinite | AbnormalFraction::NaN => {
                false
            }
        }
    }
}

impl One for AbnormalFraction {
    fn one() -> Self {
        Self::Normal(Fraction::one())
    }

    fn is_one(&self) -> bool {
        match self {
            AbnormalFraction::Normal(f) => f.is_one(),
            AbnormalFraction::Infinite | AbnormalFraction::NegInfinite | AbnormalFraction::NaN => {
                false
            }
        }
    }
}

impl MaybeExact for AbnormalFraction {
    type Approximate = f64;
    type Exact = Rational;

    fn is_exact(&self) -> bool {
        match self {
            AbnormalFraction::Normal(f) => f.is_exact(),
            AbnormalFraction::Infinite => true,
            AbnormalFraction::NegInfinite => true,
            AbnormalFraction::NaN => true,
        }
    }

    fn extract_approx(&self) -> anyhow::Result<&Self::Approximate> {
        match self {
            AbnormalFraction::Normal(f) => f.extract_approx(),
            AbnormalFraction::Infinite | AbnormalFraction::NegInfinite | AbnormalFraction::NaN => {
                Err(anyhow!("cannot extract an approximate value"))
            }
        }
    }

    fn extract_exact(&self) -> anyhow::Result<&Self::Exact> {
        match self {
            AbnormalFraction::Normal(f) => f.extract_exact(),
            AbnormalFraction::Infinite | AbnormalFraction::NegInfinite | AbnormalFraction::NaN => {
                Err(anyhow!("cannot extract an exact value"))
            }
        }
    }
    fn to_approx(self) -> anyhow::Result<Self::Approximate> {
        match self {
            AbnormalFraction::Normal(f) => f.to_approx(),
            AbnormalFraction::Infinite | AbnormalFraction::NegInfinite | AbnormalFraction::NaN => {
                Err(anyhow!("cannot extract an approximate value"))
            }
        }
    }

    fn to_exact(self) -> anyhow::Result<Self::Exact> {
        match self {
            AbnormalFraction::Normal(f) => f.to_exact(),
            AbnormalFraction::Infinite | AbnormalFraction::NegInfinite | AbnormalFraction::NaN => {
                Err(anyhow!("cannot extract an exact value"))
            }
        }
    }
}

impl PartialOrd for AbnormalFraction {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match (self, other) {
            (AbnormalFraction::Normal(f1), AbnormalFraction::Normal(f2)) => f1.partial_cmp(f2),
            (AbnormalFraction::Normal(_), AbnormalFraction::Infinite) => Some(Ordering::Less),
            (AbnormalFraction::Normal(_), AbnormalFraction::NegInfinite) => Some(Ordering::Greater),
            (AbnormalFraction::Infinite, AbnormalFraction::Normal(_)) => Some(Ordering::Greater),
            (AbnormalFraction::Infinite, AbnormalFraction::Infinite) => None,
            (AbnormalFraction::Infinite, AbnormalFraction::NegInfinite) => Some(Ordering::Greater),
            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(_)) => Some(Ordering::Less),
            (AbnormalFraction::NegInfinite, AbnormalFraction::Infinite) => Some(Ordering::Less),
            (AbnormalFraction::NegInfinite, AbnormalFraction::NegInfinite) => None,
            (_, AbnormalFraction::NaN) => None,
            (AbnormalFraction::NaN, _) => None,
        }
    }
}

impl Round for AbnormalFraction {
    fn floor(self) -> Self {
        match self {
            AbnormalFraction::Normal(f) => Self::Normal(f.floor()),
            AbnormalFraction::Infinite => Self::Infinite,
            AbnormalFraction::NegInfinite => Self::NegInfinite,
            AbnormalFraction::NaN => Self::NaN,
        }
    }

    fn ceil(self) -> Self {
        match self {
            AbnormalFraction::Normal(f) => Self::Normal(f.ceil()),
            AbnormalFraction::Infinite => Self::Infinite,
            AbnormalFraction::NegInfinite => Self::NegInfinite,
            AbnormalFraction::NaN => Self::NaN,
        }
    }
}

impl Neg for AbnormalFraction {
    type Output = AbnormalFraction;

    fn neg(self) -> Self::Output {
        match self {
            AbnormalFraction::Normal(f) => Self::Normal(-f),
            AbnormalFraction::Infinite => Self::NegInfinite,
            AbnormalFraction::NegInfinite => Self::Infinite,
            AbnormalFraction::NaN => Self::NaN,
        }
    }
}

impl Neg for &AbnormalFraction {
    type Output = AbnormalFraction;

    fn neg(self) -> Self::Output {
        match self {
            AbnormalFraction::Normal(f) => AbnormalFraction::Normal(-f),
            AbnormalFraction::Infinite => AbnormalFraction::NegInfinite,
            AbnormalFraction::NegInfinite => AbnormalFraction::Infinite,
            AbnormalFraction::NaN => AbnormalFraction::NaN,
        }
    }
}

impl Add for AbnormalFraction {
    type Output = AbnormalFraction;

    fn add(self, rhs: Self) -> Self::Output {
        print!("add {} + {}", self, rhs);
        let x = match (self, rhs) {
            (AbnormalFraction::Normal(f1), AbnormalFraction::Normal(f2)) => Self::Normal(f1 + f2),
            (AbnormalFraction::Normal(_), AbnormalFraction::Infinite) => AbnormalFraction::Infinite,
            (AbnormalFraction::Normal(_), AbnormalFraction::NegInfinite) => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::Infinite, AbnormalFraction::Normal(_)) => AbnormalFraction::Infinite,
            (AbnormalFraction::Infinite, AbnormalFraction::Infinite) => AbnormalFraction::Infinite,
            (AbnormalFraction::Infinite, AbnormalFraction::NegInfinite) => AbnormalFraction::NaN,
            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(_)) => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Infinite) => AbnormalFraction::NaN,
            (AbnormalFraction::NegInfinite, AbnormalFraction::NegInfinite) => {
                AbnormalFraction::NegInfinite
            }
            (_, AbnormalFraction::NaN) => AbnormalFraction::NaN,
            (AbnormalFraction::NaN, _) => AbnormalFraction::NaN,
        };
        println!(" = {}", x);
        x
    }
}

impl Add for &AbnormalFraction {
    type Output = AbnormalFraction;

    fn add(self, rhs: Self) -> Self::Output {
        let x = match (self, rhs) {
            (AbnormalFraction::Normal(f1), AbnormalFraction::Normal(f2)) => {
                AbnormalFraction::Normal(f1 + f2)
            }
            (AbnormalFraction::Normal(_), AbnormalFraction::Infinite) => AbnormalFraction::Infinite,
            (AbnormalFraction::Normal(_), AbnormalFraction::NegInfinite) => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::Infinite, AbnormalFraction::Normal(_)) => AbnormalFraction::Infinite,
            (AbnormalFraction::Infinite, AbnormalFraction::Infinite) => AbnormalFraction::Infinite,
            (AbnormalFraction::Infinite, AbnormalFraction::NegInfinite) => AbnormalFraction::NaN,
            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(_)) => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Infinite) => AbnormalFraction::NaN,
            (AbnormalFraction::NegInfinite, AbnormalFraction::NegInfinite) => {
                AbnormalFraction::NegInfinite
            }
            (_, AbnormalFraction::NaN) => AbnormalFraction::NaN,
            (AbnormalFraction::NaN, _) => AbnormalFraction::NaN,
        };
        println!("add {} + {} = {}", self, rhs, x);
        x
    }
}

impl AddAssign for AbnormalFraction {
    fn add_assign(&mut self, rhs: Self) {
        print!("add_assign {} + {}", self, rhs);
        if self.both_normal(&rhs) {
            if let (AbnormalFraction::Normal(f1), AbnormalFraction::Normal(f2)) = (self, rhs) {
                *f1 += f2;

                println!(" = {}", f1);
            } else {
                unreachable!()
            }
        } else {
            match (&self, &rhs) {
                (AbnormalFraction::Normal(_), AbnormalFraction::Normal(_)) => unreachable!(),
                (AbnormalFraction::Normal(_), AbnormalFraction::Infinite) => {
                    *self = AbnormalFraction::Infinite;
                }
                (AbnormalFraction::Normal(_), AbnormalFraction::NegInfinite) => {
                    *self = AbnormalFraction::NegInfinite;
                }
                (AbnormalFraction::Infinite, AbnormalFraction::Normal(_)) => {}
                (AbnormalFraction::Infinite, AbnormalFraction::Infinite) => {}
                (AbnormalFraction::Infinite, AbnormalFraction::NegInfinite) => {
                    *self = AbnormalFraction::NaN;
                }
                (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(_)) => {}
                (AbnormalFraction::NegInfinite, AbnormalFraction::Infinite) => {
                    *self = AbnormalFraction::NaN;
                }
                (AbnormalFraction::NegInfinite, AbnormalFraction::NegInfinite) => {}
                (_, AbnormalFraction::NaN) => *self = AbnormalFraction::NaN,
                (AbnormalFraction::NaN, _) => {}
            };
            println!(" = {}", self);
        }
    }
}

impl Sub for AbnormalFraction {
    type Output = AbnormalFraction;

    fn sub(self, rhs: Self) -> Self::Output {
        println!("sub {} - {}", self, rhs);
        match (&self, &rhs) {
            (AbnormalFraction::Normal(f1), AbnormalFraction::Normal(f2)) => {
                AbnormalFraction::Normal(f1 - f2)
            }
            (AbnormalFraction::Normal(_), AbnormalFraction::Infinite) => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::Normal(_), AbnormalFraction::NegInfinite) => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::Infinite, AbnormalFraction::Normal(_)) => AbnormalFraction::Infinite,
            (AbnormalFraction::Infinite, AbnormalFraction::Infinite) => AbnormalFraction::NaN,
            (AbnormalFraction::Infinite, AbnormalFraction::NegInfinite) => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(_)) => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Infinite) => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::NegInfinite) => AbnormalFraction::NaN,
            (_, AbnormalFraction::NaN) => AbnormalFraction::NaN,
            (AbnormalFraction::NaN, _) => AbnormalFraction::NaN,
        }
    }
}

impl Sub for &AbnormalFraction {
    type Output = AbnormalFraction;

    fn sub(self, rhs: Self) -> Self::Output {
        println!("sub {} - {}", self, rhs);
        match (&self, &rhs) {
            (AbnormalFraction::Normal(f1), AbnormalFraction::Normal(f2)) => {
                AbnormalFraction::Normal(f1 - f2)
            }
            (AbnormalFraction::Normal(_), AbnormalFraction::Infinite) => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::Normal(_), AbnormalFraction::NegInfinite) => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::Infinite, AbnormalFraction::Normal(_)) => AbnormalFraction::Infinite,
            (AbnormalFraction::Infinite, AbnormalFraction::Infinite) => AbnormalFraction::NaN,
            (AbnormalFraction::Infinite, AbnormalFraction::NegInfinite) => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(_)) => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Infinite) => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::NegInfinite) => AbnormalFraction::NaN,
            (_, AbnormalFraction::NaN) => AbnormalFraction::NaN,
            (AbnormalFraction::NaN, _) => AbnormalFraction::NaN,
        }
    }
}

impl SubAssign for AbnormalFraction {
    fn sub_assign(&mut self, rhs: Self) {
        println!("sub_assign {} - {}", self, rhs);
        if self.both_normal(&rhs) {
            if let (AbnormalFraction::Normal(f1), AbnormalFraction::Normal(f2)) = (self, rhs) {
                *f1 -= f2;
            } else {
                unreachable!()
            }
        } else {
            match (&self, &rhs) {
                (AbnormalFraction::Normal(_), AbnormalFraction::Normal(_)) => unreachable!(),
                (AbnormalFraction::Normal(_), AbnormalFraction::Infinite) => {
                    *self = AbnormalFraction::NegInfinite;
                }
                (AbnormalFraction::Normal(_), AbnormalFraction::NegInfinite) => {
                    *self = AbnormalFraction::Infinite;
                }
                (AbnormalFraction::Infinite, AbnormalFraction::Normal(_)) => {}
                (AbnormalFraction::Infinite, AbnormalFraction::Infinite) => {
                    *self = AbnormalFraction::NaN;
                }
                (AbnormalFraction::Infinite, AbnormalFraction::NegInfinite) => {}
                (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(_)) => {}
                (AbnormalFraction::NegInfinite, AbnormalFraction::Infinite) => {}
                (AbnormalFraction::NegInfinite, AbnormalFraction::NegInfinite) => {
                    *self = AbnormalFraction::NaN;
                }
                (_, AbnormalFraction::NaN) => *self = AbnormalFraction::NaN,
                (AbnormalFraction::NaN, _) => {}
            };
        }
    }
}

impl Mul for AbnormalFraction {
    type Output = AbnormalFraction;

    fn mul(self, rhs: Self) -> Self::Output {
        print!("mul {} * {}", self, rhs);
        let x = match (&self, &rhs) {
            (AbnormalFraction::Normal(f1), AbnormalFraction::Normal(f2)) => {
                AbnormalFraction::Normal(f1 * f2)
            }
            (AbnormalFraction::Normal(f), AbnormalFraction::Infinite) if f.is_positive() => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::Normal(f), AbnormalFraction::Infinite) if f.is_negative() => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::Normal(_), AbnormalFraction::Infinite) => AbnormalFraction::NaN,

            (AbnormalFraction::Normal(f), AbnormalFraction::NegInfinite) if f.is_positive() => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::Normal(f), AbnormalFraction::NegInfinite) if f.is_negative() => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::Normal(_), AbnormalFraction::NegInfinite) => AbnormalFraction::NaN,

            (AbnormalFraction::Infinite, AbnormalFraction::Normal(f)) if f.is_positive() => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::Infinite, AbnormalFraction::Normal(f)) if f.is_negative() => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::Infinite, AbnormalFraction::Normal(_)) => AbnormalFraction::NaN,

            (AbnormalFraction::Infinite, AbnormalFraction::Infinite) => AbnormalFraction::Infinite,
            (AbnormalFraction::Infinite, AbnormalFraction::NegInfinite) => {
                AbnormalFraction::NegInfinite
            }

            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(f)) if f.is_positive() => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(f)) if f.is_negative() => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(_)) => AbnormalFraction::NaN,

            (AbnormalFraction::NegInfinite, AbnormalFraction::Infinite) => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::NegInfinite) => {
                AbnormalFraction::Infinite
            }
            (_, AbnormalFraction::NaN) => AbnormalFraction::NaN,
            (AbnormalFraction::NaN, _) => AbnormalFraction::NaN,
        };
        println!(" = {}", x);
        x
    }
}

impl Mul for &AbnormalFraction {
    type Output = AbnormalFraction;

    fn mul(self, rhs: Self) -> Self::Output {
        print!("mul {} * {}", self, rhs);
        let x = match (&self, &rhs) {
            (AbnormalFraction::Normal(f1), AbnormalFraction::Normal(f2)) => {
                AbnormalFraction::Normal(f1 * f2)
            }
            (AbnormalFraction::Normal(f), AbnormalFraction::Infinite) if f.is_positive() => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::Normal(f), AbnormalFraction::Infinite) if f.is_negative() => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::Normal(_), AbnormalFraction::Infinite) => AbnormalFraction::NaN,

            (AbnormalFraction::Normal(f), AbnormalFraction::NegInfinite) if f.is_positive() => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::Normal(f), AbnormalFraction::NegInfinite) if f.is_negative() => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::Normal(_), AbnormalFraction::NegInfinite) => AbnormalFraction::NaN,

            (AbnormalFraction::Infinite, AbnormalFraction::Normal(f)) if f.is_positive() => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::Infinite, AbnormalFraction::Normal(f)) if f.is_negative() => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::Infinite, AbnormalFraction::Normal(_)) => AbnormalFraction::NaN,

            (AbnormalFraction::Infinite, AbnormalFraction::Infinite) => AbnormalFraction::Infinite,
            (AbnormalFraction::Infinite, AbnormalFraction::NegInfinite) => {
                AbnormalFraction::NegInfinite
            }

            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(f)) if f.is_positive() => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(f)) if f.is_negative() => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(_)) => AbnormalFraction::NaN,

            (AbnormalFraction::NegInfinite, AbnormalFraction::Infinite) => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::NegInfinite) => {
                AbnormalFraction::Infinite
            }
            (_, AbnormalFraction::NaN) => AbnormalFraction::NaN,
            (AbnormalFraction::NaN, _) => AbnormalFraction::NaN,
        };
        println!(" = {}", x);
        x
    }
}

impl Div for AbnormalFraction {
    type Output = AbnormalFraction;

    fn div(self, rhs: Self) -> Self::Output {
        print!("div {} / {}", self, rhs);
        let x = match (&self, &rhs) {
            (AbnormalFraction::Normal(f1), AbnormalFraction::Normal(f2)) if !f2.is_zero() => {
                AbnormalFraction::Normal(f1 / f2)
            }
            (AbnormalFraction::Normal(_), AbnormalFraction::Normal(_)) => AbnormalFraction::NaN,
            (AbnormalFraction::Normal(_), AbnormalFraction::Infinite) => {
                AbnormalFraction::Normal(Fraction::zero())
            }
            (AbnormalFraction::Normal(_), AbnormalFraction::NegInfinite) => {
                AbnormalFraction::Normal(Fraction::zero())
            }

            (AbnormalFraction::Infinite, AbnormalFraction::Normal(f)) if f.is_positive() => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::Infinite, AbnormalFraction::Normal(f)) if f.is_negative() => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::Infinite, AbnormalFraction::Normal(_)) => AbnormalFraction::NaN,

            (AbnormalFraction::Infinite, AbnormalFraction::Infinite) => AbnormalFraction::NaN,
            (AbnormalFraction::Infinite, AbnormalFraction::NegInfinite) => AbnormalFraction::NaN,

            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(f)) if f.is_positive() => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(f)) if f.is_negative() => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(_)) => AbnormalFraction::NaN,

            (AbnormalFraction::NegInfinite, AbnormalFraction::Infinite) => AbnormalFraction::NaN,
            (AbnormalFraction::NegInfinite, AbnormalFraction::NegInfinite) => AbnormalFraction::NaN,
            (_, AbnormalFraction::NaN) => AbnormalFraction::NaN,
            (AbnormalFraction::NaN, _) => AbnormalFraction::NaN,
        };
        println!(" = {}", x);
        x
    }
}

impl Div for &AbnormalFraction {
    type Output = AbnormalFraction;

    fn div(self, rhs: Self) -> Self::Output {
        print!("div {} / {}", self, rhs);
        let x = match (&self, &rhs) {
            (AbnormalFraction::Normal(f1), AbnormalFraction::Normal(f2)) if !f2.is_zero() => {
                AbnormalFraction::Normal(f1 / f2)
            }
            (AbnormalFraction::Normal(_), AbnormalFraction::Normal(_)) => AbnormalFraction::NaN,
            (AbnormalFraction::Normal(_), AbnormalFraction::Infinite) => {
                AbnormalFraction::Normal(Fraction::zero())
            }
            (AbnormalFraction::Normal(_), AbnormalFraction::NegInfinite) => {
                AbnormalFraction::Normal(Fraction::zero())
            }

            (AbnormalFraction::Infinite, AbnormalFraction::Normal(f)) if f.is_positive() => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::Infinite, AbnormalFraction::Normal(f)) if f.is_negative() => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::Infinite, AbnormalFraction::Normal(_)) => AbnormalFraction::NaN,

            (AbnormalFraction::Infinite, AbnormalFraction::Infinite) => AbnormalFraction::NaN,
            (AbnormalFraction::Infinite, AbnormalFraction::NegInfinite) => AbnormalFraction::NaN,

            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(f)) if f.is_positive() => {
                AbnormalFraction::NegInfinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(f)) if f.is_negative() => {
                AbnormalFraction::Infinite
            }
            (AbnormalFraction::NegInfinite, AbnormalFraction::Normal(_)) => AbnormalFraction::NaN,

            (AbnormalFraction::NegInfinite, AbnormalFraction::Infinite) => AbnormalFraction::NaN,
            (AbnormalFraction::NegInfinite, AbnormalFraction::NegInfinite) => AbnormalFraction::NaN,
            (_, AbnormalFraction::NaN) => AbnormalFraction::NaN,
            (AbnormalFraction::NaN, _) => AbnormalFraction::NaN,
        };
        println!(" = {}", x);
        x
    }
}

impl Sum for AbnormalFraction {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |a, b| a + b)
    }
}

impl From<usize> for AbnormalFraction {
    fn from(value: usize) -> Self {
        Self::Normal(value.into())
    }
}

impl From<(usize, usize)> for AbnormalFraction {
    fn from(value: (usize, usize)) -> Self {
        Self::Normal(value.into())
    }
}

#[macro_export]
/// Convenience short-hand macro to create fractions.
macro_rules! f_ab {
    ($e: expr) => {
        AbnormalFraction::from($e)
    };

    ($e: expr, $f: expr) => {
        AbnormalFraction::from(($e, $f))
    };
}
pub use f_ab;

#[macro_export]
/// Convenience short-hand macro to create a fraction representing zero.
macro_rules! f0_ab {
    () => {
        AbnormalFraction::zero()
    };
}
pub use f0_ab;

#[macro_export]
/// Convenience short-hand macro to create a fraction representing one.
macro_rules! f1_ab {
    () => {
        AbnormalFraction::one()
    };
}
pub use f1_ab;
use malachite::rational::Rational;
use pathfinding::num_traits;

#[cfg(test)]
mod tests {
    use ebi_arithmetic::Zero;

    use crate::abnormal_fraction::AbnormalFraction;

    #[test]
    fn abnormal_fraction() {
        assert!(AbnormalFraction::zero().is_zero());
        assert!(!AbnormalFraction::infinity().is_zero());
        assert!(AbnormalFraction::infinity().is_infinite());
        assert!(!AbnormalFraction::infinity().is_finite());
    }
}
