use crate::raytrace::vec3::{Point3, Vec3};

#[derive(Debug, Copy, Clone, PartialEq, Default)]
pub struct Ray {
    pub orig: Point3,
    pub dir: Vec3,
    pub time: f64,
}

impl Ray {
    pub fn new(orig: Point3, dir: Vec3, time: f64) -> Ray {
        Ray { orig, dir, time }
    }

    pub fn at(&self, t: f64) -> Point3 {
        self.orig + t * self.dir
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn basic() {
        let orig = Point3::new(3.0, 3.0, 4.0);
        let dir = Vec3::new(2.0, 6.0, 4.0);
        let r = Ray::new(orig, dir);

        assert_eq!(r.at(0.0), orig);
        assert_eq!(r.at(1.0), orig + dir);
    }
}
