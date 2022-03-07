use super::{
    ray::Ray,
    vec3::{cross, unit, Point3, Vec3},
};

pub struct Camera {
    pub origin: Point3,
    pub lower_left_corner: Point3,
    pub horizontal: Vec3,
    pub vertical: Vec3,
}

impl Camera {
    pub fn new(lookfrom: Point3, lookat: Point3, vup: Vec3, fov: f64, aspect_ratio: f64) -> Camera {
        let theta = fov.to_radians();
        let h = f64::tan(theta / 2.0);
        let viewport_height = 2.0 * h;
        let viewport_width = aspect_ratio * viewport_height;

        let w = unit(lookfrom - lookat);
        let u = unit(cross(vup, w));
        let v = cross(w, u);

        let origin = lookfrom;
        let horizontal = viewport_width * u;
        let vertical = viewport_height * v;
        let lower_left_corner = origin - horizontal / 2.0 - vertical / 2.0 - w;

        Camera {
            origin,
            lower_left_corner,
            horizontal,
            vertical,
        }
    }

    pub fn get_ray(&self, u: f64, v: f64) -> Ray {
        let origin = self.origin;
        let direction =
            self.lower_left_corner + u * self.horizontal + v * self.vertical - self.origin;

        Ray::new(origin, direction)
    }
}
