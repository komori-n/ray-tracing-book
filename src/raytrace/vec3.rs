use rand::distributions::Distribution;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Debug, Copy, Clone, PartialEq, Default)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
pub type Point3 = Vec3;

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3 { x, y, z }
    }

    pub fn length(self) -> f64 {
        f64::sqrt(self.length_squared())
    }

    pub fn length_squared(self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn random(rng: &mut dyn rand::RngCore) -> Vec3 {
        Self::random_with_minmax(rng, 0.0, 1.0)
    }

    pub fn random_with_minmax(rng: &mut dyn rand::RngCore, min: f64, max: f64) -> Vec3 {
        let uni = rand::distributions::Uniform::from(min..max);
        let x = uni.sample(rng);
        let y = uni.sample(rng);
        let z = uni.sample(rng);
        Vec3::new(x, y, z)
    }
}

impl std::fmt::Display for Vec3 {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} {} {}", self.x, self.y, self.z)
    }
}

pub fn dot(lhs: Vec3, rhs: Vec3) -> f64 {
    lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z
}

pub fn cross(lhs: Vec3, rhs: Vec3) -> Vec3 {
    Vec3 {
        x: lhs.y * rhs.z - lhs.z * rhs.y,
        y: lhs.z * rhs.x - lhs.x * rhs.z,
        z: lhs.x * rhs.y - lhs.y * rhs.x,
    }
}

pub fn unit(v: Vec3) -> Vec3 {
    let len = v.length();
    Vec3 {
        x: v.x / len,
        y: v.y / len,
        z: v.z / len,
    }
}

pub fn random_in_unit_sphere(rng: &mut dyn rand::RngCore) -> Vec3 {
    loop {
        let p = Vec3::random_with_minmax(rng, -1.0, 1.0);
        if p.length_squared() < 1.0 {
            return p;
        }
    }
}

pub fn random_unit_veector(rng: &mut dyn rand::RngCore) -> Vec3 {
    unit(random_in_unit_sphere(rng))
}

pub fn random_in_hemisphere(rng: &mut dyn rand::RngCore, n: &Vec3) -> Vec3 {
    let in_unit_sphere = random_in_unit_sphere(rng);
    if dot(in_unit_sphere, *n) > 0.0 {
        in_unit_sphere
    } else {
        -in_unit_sphere
    }
}

impl Add for Vec3 {
    type Output = Vec3;

    fn add(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for Vec3 {
    type Output = Vec3;

    fn sub(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Mul<Vec3> for Vec3 {
    type Output = f64;

    fn mul(self, other: Vec3) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;

    fn mul(self, other: f64) -> Vec3 {
        Vec3 {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;

    fn mul(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self * other.x,
            y: self * other.y,
            z: self * other.z,
        }
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;

    fn div(self, other: f64) -> Vec3 {
        Vec3 {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

impl Neg for Vec3 {
    type Output = Vec3;

    fn neg(self) -> Vec3 {
        Vec3 {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, other: Vec3) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl SubAssign for Vec3 {
    fn sub_assign(&mut self, other: Vec3) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl MulAssign<f64> for Vec3 {
    fn mul_assign(&mut self, other: f64) {
        self.x *= other;
        self.y *= other;
        self.z *= other;
    }
}

impl DivAssign<f64> for Vec3 {
    fn div_assign(&mut self, other: f64) {
        self.x /= other;
        self.y /= other;
        self.z /= other;
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn basic() {
        let a = Vec3::new(3.0, 3.0, 4.0);
        let b = Vec3::new(2.0, 6.0, 4.0);

        assert_eq!(a.to_string(), "3 3 4".to_string());
        assert_eq!(a.length_squared(), 3.0 * 3.0 + 3.0 * 3.0 + 4.0 * 4.0);
        assert_eq!(a.length(), f64::sqrt(3.0 * 3.0 + 3.0 * 3.0 + 4.0 * 4.0));

        assert_eq!(dot(a, b), 40.0);
        assert_eq!(cross(a, b), Vec3::new(-12.0, -4.0, 12.0));
        assert_eq!(
            unit(a),
            Vec3::new(3.0, 3.0, 4.0) / f64::sqrt(3.0 * 3.0 + 3.0 * 3.0 + 4.0 * 4.0)
        );
        assert_eq!(a + b, Vec3::new(5.0, 9.0, 8.0));
        assert_eq!(a - b, Vec3::new(1.0, -3.0, 0.0));
        assert_eq!(a * b, 40.0);
        assert_eq!(a * 2.0, Vec3::new(6.0, 6.0, 8.0));
        assert_eq!(2.0 * a, Vec3::new(6.0, 6.0, 8.0));
        assert_eq!(b / 2.0, Vec3::new(1.0, 3.0, 2.0));
        assert_eq!(-a, Vec3::new(-3.0, -3.0, -4.0));

        let mut c = a;
        c += b;
        assert_eq!(c, Vec3::new(5.0, 9.0, 8.0));

        let mut c = a;
        c -= b;
        assert_eq!(c, Vec3::new(1.0, -3.0, 0.0));

        let mut c = a;
        c *= 2.0;
        assert_eq!(c, Vec3::new(6.0, 6.0, 8.0));

        let mut c = b;
        c /= 2.0;
        assert_eq!(c, Vec3::new(1.0, 3.0, 2.0));
    }
}
