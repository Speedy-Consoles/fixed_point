mod precalc;

#[macro_use] extern crate derive_more;

use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Sub;
use std::ops::SubAssign;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Div;
use std::ops::DivAssign;
use std::f32::consts::PI as PI32;
use std::f64::consts::PI as PI64;
use std::fmt;

use cgmath::Vector3;
use cgmath::Vector2;
use cgmath::Point3;
use cgmath::Point2;

use self::precalc::SIN_PRECISION;
use self::precalc::SIN;

// fixed-point constants
const FP_PRECISION: u64 = 16;
const FP_RESOLUTION: u64 = 1 << FP_PRECISION;

// trigonometry constants
const FP_SIN_PRECISION_DIFF: u64 = FP_PRECISION - SIN_PRECISION;
const SIN_QUARTER_RESOLUTION: u64 = 1 << (SIN_PRECISION - 2);
const FP_SIN_RESOLUTION_RATIO: u64 = 1 << FP_SIN_PRECISION_DIFF;
const SIN_QUARTER_MASK: u64 = (!0) % SIN_QUARTER_RESOLUTION;
const SIN_INTRA_MASK: u64 = (!0) % FP_SIN_RESOLUTION_RATIO;

#[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord,
    Add, AddAssign, Sub, SubAssign, Rem, RemAssign, Neg)]
pub struct FixedPoint(i64);


// number of bits in the fractional part: FP_PRECISION
// number of bits in the whole part: 64 - 2 * FP_PRECISION
impl FixedPoint {
    pub const fn new(value: i64) -> FixedPoint {
        FixedPoint(value << FP_PRECISION)
    }

    pub const fn zero() -> FixedPoint {
        FixedPoint(0)
    }

    pub const fn one() -> FixedPoint {
        FixedPoint(FP_RESOLUTION as i64)
    }

    pub const fn fraction(nominator: i64, denominator: i64) -> Self {
        FixedPoint((nominator << FP_PRECISION) / denominator)
    }

    // TODO make const
    pub fn abs(&self) -> FixedPoint {
        FixedPoint(self.0.abs())
    }

    // TODO make const
    pub fn inv_sqrt(&self) -> FixedPoint {
        const THREE: i64 = FP_RESOLUTION as i64 * 3;
        if self.0 <= 0 {
            panic!("Attempted to take inverse square root of non-positive number!");
        }
        // TODO find better approximation
        //eprintln!("starting inverse square root computation of {:?}", self);
        let shift = 64 - self.0.leading_zeros() as i64 - FP_PRECISION as i64;
        //eprintln!("shift: {}", shift);
        //eprintln!("sqrt_approx: {} ({})", FixedPoint(FP_PRECISION as i64 + shift / 2), Into::<f64>::into(*self).sqrt());
        let mut approx = (1 << FP_PRECISION as i64 - shift / 2);
        //eprintln!("approx: {} ({})", FixedPoint(approx), 1.0 / Into::<f64>::into(*self).sqrt());
        for _ in 0..5 { // TODO relate number of iterations to FP_PRECISION
            approx = fp_mul(THREE - fp_mul(fp_mul(self.0, approx), approx), approx) >> 1;
        }
        FixedPoint(approx)
    }

    pub const fn is_zero(&self) -> bool {
        self.0 == 0
    }

    pub const fn is_positive(&self) -> bool {
        self.0 > 0
    }

    pub const fn is_negative(&self) -> bool {
        self.0 < 0
    }

    pub const fn mix(self, other: FixedPoint, ratio: FixedPoint) -> FixedPoint {
        FixedPoint(fp_mul(self.0, FP_RESOLUTION as i64 - ratio.0) + fp_mul(other.0, ratio.0))
    }
}

const fn fp_mul(a: i64, b: i64) -> i64 {
    (a * b) >> FP_PRECISION
}

const fn fp_div(a: i64, b: i64) -> i64 {
    (a << FP_PRECISION) / b
}

impl Mul<FixedPoint> for FixedPoint {
    type Output = FixedPoint;

    fn mul(self, rhs: FixedPoint) -> FixedPoint {
        FixedPoint(fp_mul(self.0, rhs.0))
    }
}

impl Mul<i64> for FixedPoint {
    type Output = FixedPoint;

    fn mul(self, rhs: i64) -> FixedPoint {
        FixedPoint(self.0 * rhs)
    }
}

impl MulAssign<FixedPoint> for FixedPoint {
    fn mul_assign(&mut self, rhs: FixedPoint) {
        self.0 = fp_mul(self.0, rhs.0);
    }
}

impl MulAssign<i64> for FixedPoint {
    fn mul_assign(&mut self, rhs: i64) {
        self.0 *= rhs;
    }
}

impl Div<FixedPoint> for FixedPoint {
    type Output = FixedPoint;

    fn div(self, rhs: FixedPoint) -> FixedPoint {
        FixedPoint(fp_div(self.0, rhs.0))
    }
}

impl Div<i64> for FixedPoint {
    type Output = FixedPoint;

    fn div(self, rhs: i64) -> FixedPoint {
        FixedPoint(self.0 / rhs)
    }
}

impl DivAssign<FixedPoint> for FixedPoint {
    fn div_assign(&mut self, rhs: FixedPoint) {
        self.0 = fp_div(self.0, rhs.0);
    }
}

impl DivAssign<i64> for FixedPoint {
    fn div_assign(&mut self, rhs: i64) {
        self.0 /= rhs;
    }
}

impl From<i64> for FixedPoint {
    fn from(value: i64) -> Self {
        FixedPoint::new(value)
    }
}

impl From<f64> for FixedPoint {
    fn from(value: f64) -> Self {
        FixedPoint((value * FP_RESOLUTION as f64).round() as i64)
    }
}

impl From<f32> for FixedPoint {
    fn from(value: f32) -> Self {
        FixedPoint((value as f64 * FP_RESOLUTION as f64).round() as i64)
    }
}

impl Into<f64> for FixedPoint {
    fn into(self) -> f64 {
        self.0 as f64 / FP_RESOLUTION as f64
    }
}

impl Into<f32> for FixedPoint {
    fn into(self) -> f32 {
        self.0 as f32 / FP_RESOLUTION as f32
    }
}

impl fmt::Display for FixedPoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let float: f64 = (*self).into();
        write!(f, "s{:.*}", (0.4 * FP_PRECISION as f64) as usize, float)
    }
}

impl fmt::Debug for FixedPoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

#[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord,
    Add, AddAssign, Sub, SubAssign, Neg, From)]
pub struct FPAngle(FixedPoint);

impl FPAngle {
    pub const fn fraction(nominator: i64, denominator: i64) -> FPAngle {
        FPAngle(FixedPoint::fraction(nominator, denominator))
    }

    pub const fn whole() -> FPAngle {
        FPAngle(FixedPoint::new(1))
    }

    pub const fn half() -> FPAngle {
        FPAngle(FixedPoint::fraction(1, 2))
    }

    pub const fn quarter() -> FPAngle {
        FPAngle(FixedPoint::fraction(1, 4))
    }

    pub const fn zero() -> FPAngle {
        FPAngle(FixedPoint::new(0))
    }

    // TODO make const
    pub fn sin(&self) -> FixedPoint {
        let circular = (((self.0).0 % FP_RESOLUTION as i64)
            + FP_RESOLUTION as i64) as u64 % FP_RESOLUTION;
        let intra = circular & SIN_INTRA_MASK;
        let quadrant = circular >> (FP_PRECISION - 2);
        let mut index = ((circular >> FP_SIN_PRECISION_DIFF) & SIN_QUARTER_MASK) as usize;
        let mut next_index = index + 1;
        if quadrant % 2 != 0 {
            index = SIN_QUARTER_RESOLUTION as usize - index;
            next_index = index - 1;
        };
        let sin1 = SIN[index];
        let sin2 = SIN[next_index];
        let mut sin = (sin1 * (FP_SIN_RESOLUTION_RATIO - intra) as i64
            + sin2 * intra as i64) >> FP_SIN_PRECISION_DIFF as i64;
        if quadrant / 2 != 0 {
            sin = -sin
        };
        FixedPoint(sin)
    }

    // TODO make const
    pub fn cos(&self) -> FixedPoint {
        (*self + FPAngle::quarter()).sin()
    }

    pub fn from_tau_float(float: f64) -> FPAngle {
        FPAngle(FixedPoint((float * FP_RESOLUTION as f64) as i64))
    }

    pub fn rad_f32(self) -> f32 {
        let f: f32 = self.0.into();
        f * 2.0 * PI32
    }

    pub fn rad_f64(self) -> f64 {
        let f: f64 = self.0.into();
        f * 2.0 * PI64
    }
}

impl Mul<FixedPoint> for FPAngle {
    type Output = Self;

    fn mul(self, rhs: FixedPoint) -> Self {
        FPAngle(self.0 * rhs.0)
    }
}

impl fmt::Display for FPAngle {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let float: f64 = self.0.into();
        write!(f, "a{:.*}", (0.3 * FP_PRECISION as f64) as usize, float)
    }
}

impl fmt::Debug for FPAngle {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

#[derive(Clone, Copy)]
pub struct FPVec2 {
    pub x: FixedPoint,
    pub y: FixedPoint,
}

impl FPVec2 {
    pub const fn new(x: FixedPoint, y: FixedPoint) -> FPVec2 {
        FPVec2 { x, y }
    }

    pub const fn zero() -> FPVec2 {
        FPVec2 { x: FixedPoint::zero(), y: FixedPoint::zero() }
    }

    pub const fn dot(&self, other: FPVec2) -> FixedPoint {
        FixedPoint(fp_mul(self.x.0, other.x.0) + fp_mul(self.y.0, other.y.0))
    }

    pub const fn length2(&self) -> FixedPoint {
        self.dot(*self)
    }

    // TODO make const
    pub fn scale_to(&self, length: FixedPoint) -> FPVec2 {
        let fp_length = fp_mul(length.0, self.length2().inv_sqrt().0);
        FPVec2 {
            x: FixedPoint(self.x.0 * fp_length),
            y: FixedPoint(self.y.0 * fp_length),
        }
    }

    // TODO make const
    pub fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero()
    }
}

impl Add for FPVec2 {
    type Output = FPVec2;

    fn add(self, other: FPVec2) -> Self::Output {
        FPVec2 {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Sub for FPVec2 {
    type Output = FPVec2;

    fn sub(self, other: Self) -> FPVec2 {
        FPVec2 {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Mul<FixedPoint> for FPVec2 {
    type Output = FPVec2;

    fn mul(self, rhs: FixedPoint) -> FPVec2 {
        FPVec2 {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl Mul<FPVec2> for FixedPoint {
    type Output = FPVec2;

    fn mul(self, rhs: FPVec2) -> FPVec2 {
        FPVec2 {
            x: self * rhs.x,
            y: self * rhs.y,
        }
    }
}

impl MulAssign<FixedPoint> for FPVec2 {
    fn mul_assign(&mut self, rhs: FixedPoint) {
        self.x *= rhs;
        self.y *= rhs;
    }
}

impl Div<FixedPoint> for FPVec2 {
    type Output = FPVec2;

    fn div(self, rhs: FixedPoint) -> FPVec2 {
        FPVec2 {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

impl DivAssign<FixedPoint> for FPVec2 {
    fn div_assign(&mut self, rhs: FixedPoint) {
        self.x /= rhs;
        self.y /= rhs;
    }
}

impl AddAssign for FPVec2 {
    fn add_assign(&mut self, other: FPVec2) {
        self.x += other.x;
        self.y += other.y;
    }
}

impl SubAssign for FPVec2 {
    fn sub_assign(&mut self, other: FPVec2) {
        self.x -= other.x;
        self.y -= other.y;
    }
}

impl Into<Vector2<f32>> for FPVec2 {
    fn into(self) -> Vector2<f32> {
        Vector2::new(
            self.x.into(),
            self.y.into(),
        )
    }
}

impl Into<Point2<f32>> for FPVec2 {
    fn into(self) -> Point2<f32> {
        Point2::new(
            self.x.into(),
            self.y.into(),
        )
    }
}

impl fmt::Display for FPVec2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "v({}, {})", self.x, self.y)
    }
}

impl fmt::Debug for FPVec2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

#[derive(Clone, Copy)]
pub struct FPVec3 {
    pub x: FixedPoint,
    pub y: FixedPoint,
    pub z: FixedPoint,
}

impl FPVec3 {
    pub const fn new(x: FixedPoint, y: FixedPoint, z: FixedPoint) -> FPVec3 {
        FPVec3 { x, y, z }
    }

    pub const fn zero() -> FPVec3 {
        FPVec3 { x: FixedPoint::zero(), y: FixedPoint::zero(), z: FixedPoint::zero() }
    }

    pub const fn dot(&self, other: FPVec3) -> FixedPoint {
        FixedPoint(
              fp_mul(self.x.0, other.x.0)
            + fp_mul(self.y.0, other.y.0)
            + fp_mul(self.z.0, other.z.0)
        )
    }

    pub const fn cross(&self, other: FPVec3) -> FPVec3 {
        FPVec3 {
            x: FixedPoint(fp_mul(self.y.0, other.z.0) - fp_mul(self.z.0, other.y.0)),
            y: FixedPoint(fp_mul(self.z.0, other.x.0) - fp_mul(self.x.0, other.z.0)),
            z: FixedPoint(fp_mul(self.x.0, other.y.0) - fp_mul(self.y.0, other.x.0)),
        }
    }

    pub const fn length2(&self) -> FixedPoint {
        self.dot(*self)
    }

    // TODO make const
    pub fn scale_to(&self, length: FixedPoint) -> FPVec3 {
        let fp_length = fp_mul(length.0, self.length2().inv_sqrt().0);
        FPVec3 {
            x: FixedPoint(self.x.0 * fp_length),
            y: FixedPoint(self.y.0 * fp_length),
            z: FixedPoint(self.z.0 * fp_length),
        }
    }

    // TODO make const
    pub fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero() && self.z.is_zero()
    }

    pub fn xy(&self) -> FPVec2 {
        FPVec2 { x: self.x, y: self.y }
    }
}

impl Add for FPVec3 {
    type Output = FPVec3;

    fn add(self, other: FPVec3) -> Self::Output {
        FPVec3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for FPVec3 {
    type Output = FPVec3;

    fn sub(self, other: Self) -> FPVec3 {
        FPVec3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Mul<FixedPoint> for FPVec3 {
    type Output = FPVec3;

    fn mul(self, rhs: FixedPoint) -> FPVec3 {
        FPVec3 {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Mul<FPVec3> for FixedPoint {
    type Output = FPVec3;

    fn mul(self, rhs: FPVec3) -> FPVec3 {
        FPVec3 {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
        }
    }
}

impl MulAssign<FixedPoint> for FPVec3 {
    fn mul_assign(&mut self, rhs: FixedPoint) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl Div<FixedPoint> for FPVec3 {
    type Output = FPVec3;

    fn div(self, rhs: FixedPoint) -> FPVec3 {
        FPVec3 {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl DivAssign<FixedPoint> for FPVec3 {
    fn div_assign(&mut self, rhs: FixedPoint) {
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
    }
}

impl AddAssign for FPVec3 {
    fn add_assign(&mut self, other: FPVec3) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl SubAssign for FPVec3 {
    fn sub_assign(&mut self, other: FPVec3) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl Into<Vector3<f32>> for FPVec3 {
    fn into(self) -> Vector3<f32> {
        Vector3::new(
            self.x.into(),
            self.y.into(),
            self.z.into(),
        )
    }
}

impl Into<Point3<f32>> for FPVec3 {
    fn into(self) -> Point3<f32> {
        Point3::new(
            self.x.into(),
            self.y.into(),
            self.z.into(),
        )
    }
}

impl fmt::Display for FPVec3 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "v({}, {}, {})", self.x, self.y, self.z)
    }
}

impl fmt::Debug for FPVec3 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

#[cfg(test)]
mod tests {
    use crate::FixedPoint;
    use crate::FP_PRECISION;
    use crate::FP_RESOLUTION;

    #[test]
    fn test_inv_sqrt() {
        for x in 500..40000 {
            let fixed = FixedPoint((x as f64).powf(3.1415) as i64);
            let fixed_is = fixed.inv_sqrt();
            let float = Into::<f64>::into(fixed);
            let float_is = 1.0 / Into::<f64>::into(fixed).sqrt();
            let diff = Into::<f64>::into(fixed_is) - float_is;
            assert!(
                diff.abs() < 2.0 / FP_RESOLUTION as f64,
                "Inverse square root error too big:\n{}.inv_sqrt(): {}\n1.0 / {}.sqrt(): {}\ndiff.: {}",
                fixed, fixed_is, float, float_is, diff
            );
        }
    }
}