use itertools::iproduct;
use rand::distributions::Distribution;
use std::borrow::BorrowMut;
use std::cmp::Ordering;
use std::f64::consts::PI;
use std::sync::{Arc, Mutex};

use crate::raytrace::ray::Ray;
use crate::raytrace::vec3::{dot, Point3, Vec3};

use super::material::Material;

#[derive(Debug, Default, Clone)]
pub struct AABB {
    pub min: Point3,
    pub max: Point3,
}

impl AABB {
    pub fn new(min: Point3, max: Point3) -> AABB {
        AABB { min, max }
    }

    pub fn hit(&self, r: &Ray, mut t_min: f64, mut t_max: f64) -> bool {
        for a in 0..3 {
            let inv_d = 1.0 / r.dir[a];
            let mut t0 = (self.min[a] - r.orig[a]) * inv_d;
            let mut t1 = (self.max[a] - r.orig[a]) * inv_d;
            if inv_d < 0.0 {
                std::mem::swap(&mut t0, &mut t1);
            }

            t_min = if t0 > t_min { t0 } else { t_min };
            t_max = if t1 < t_max { t1 } else { t_max };
            if t_max <= t_min {
                return false;
            }
        }
        true
    }
}

fn surrounding_box(box0: &AABB, box1: &AABB) -> AABB {
    let small = Point3::new(
        if box0.min.x < box1.min.x {
            box0.min.x
        } else {
            box1.min.x
        },
        if box0.min.y < box1.min.y {
            box0.min.y
        } else {
            box1.min.y
        },
        if box0.min.z < box1.min.z {
            box0.min.z
        } else {
            box1.min.z
        },
    );
    let big = Point3::new(
        if box0.max.x > box1.max.x {
            box0.max.x
        } else {
            box1.max.x
        },
        if box0.max.y > box1.max.y {
            box0.max.y
        } else {
            box1.max.y
        },
        if box0.max.z > box1.max.z {
            box0.max.z
        } else {
            box1.max.z
        },
    );
    AABB::new(small, big)
}

pub struct HitRecord {
    pub p: Point3,
    pub normal: Vec3,
    pub material: Arc<dyn Material + Send + Sync>,
    pub t: f64,
    pub u: f64,
    pub v: f64,
    pub front_face: bool,
}

impl HitRecord {
    fn new(
        r: &Ray,
        p: Point3,
        outward_normal: Vec3,
        material: Arc<dyn Material + Send + Sync>,
        t: f64,
        u: f64,
        v: f64,
    ) -> HitRecord {
        let front_face = dot(r.dir, outward_normal) < 0.0;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };

        HitRecord {
            p,
            normal,
            material,
            t,
            u,
            v,
            front_face,
        }
    }
}

pub trait Hittable {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB>;
}

pub struct Sphere {
    center: Point3,
    radius: f64,
    material: Arc<dyn Material + Send + Sync>,
}

impl Sphere {
    pub fn new(center: Point3, radius: f64, material: Arc<dyn Material + Send + Sync>) -> Sphere {
        Sphere {
            center,
            radius,
            material,
        }
    }

    pub fn get_sphere_uv(p: Point3) -> (f64, f64) {
        let theta = f64::acos(-p.y);
        let phi = f64::atan2(-p.z, p.x) + PI;

        (phi / (2.0 * PI), theta / PI)
    }
}

impl Hittable for Sphere {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = r.orig - self.center;
        let a = r.dir.length_squared();
        let half_b = dot(oc, r.dir);
        let c = oc.length_squared() - self.radius * self.radius;

        let discriminant = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return None;
        }

        let sqrtd = f64::sqrt(discriminant);
        let mut root = (-half_b - sqrtd) / a;
        if root < t_min || t_max < root {
            root = (-half_b + sqrtd) / a;
            if root < t_min || t_max < root {
                return None;
            }
        }
        let root = root;
        let outward_normal = (r.at(root) - self.center) / self.radius;
        let (u, v) = Sphere::get_sphere_uv(outward_normal);

        Some(HitRecord::new(
            r,
            r.at(root),
            outward_normal,
            self.material.clone(),
            root,
            u,
            v,
        ))
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        Some(AABB::new(
            self.center - Vec3::new(self.radius, self.radius, self.radius),
            self.center + Vec3::new(self.radius, self.radius, self.radius),
        ))
    }
}

pub struct MovingSphere {
    center0: Point3,
    center1: Point3,
    t0: f64,
    t1: f64,
    radius: f64,
    material: Arc<dyn Material + Send + Sync>,
}

impl MovingSphere {
    pub fn new(
        center0: Point3,
        center1: Point3,
        t0: f64,
        t1: f64,
        radius: f64,
        material: Arc<dyn Material + Send + Sync>,
    ) -> MovingSphere {
        MovingSphere {
            center0,
            center1,
            t0,
            t1,
            radius,
            material,
        }
    }

    pub fn center(&self, time: f64) -> Point3 {
        self.center0 + ((time - self.t0) / (self.t1 - self.t0)) * (self.center1 - self.center0)
    }
}

impl Hittable for MovingSphere {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = r.orig - self.center(r.time);
        let a = r.dir.length_squared();
        let half_b = dot(oc, r.dir);
        let c = oc.length_squared() - self.radius * self.radius;

        let discriminant = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return None;
        }

        let sqrtd = f64::sqrt(discriminant);
        let mut root = (-half_b - sqrtd) / a;
        if root < t_min || t_max < root {
            root = (-half_b + sqrtd) / a;
            if root < t_min || t_max < root {
                return None;
            }
        }
        let root = root;

        Some(HitRecord::new(
            r,
            r.at(root),
            (r.at(root) - self.center(r.time)) / self.radius,
            self.material.clone(),
            root,
            0.0,
            0.0,
        ))
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        let box0 = AABB::new(
            self.center(t0) - Vec3::new(self.radius, self.radius, self.radius),
            self.center(t0) + Vec3::new(self.radius, self.radius, self.radius),
        );
        let box1 = AABB::new(
            self.center(t1) - Vec3::new(self.radius, self.radius, self.radius),
            self.center(t1) + Vec3::new(self.radius, self.radius, self.radius),
        );

        Some(surrounding_box(&box0, &box1))
    }
}

pub struct BVHNode {
    left: Arc<dyn Hittable + Send + Sync>,
    right: Arc<dyn Hittable + Send + Sync>,
    box_: AABB,
}

impl BVHNode {
    pub fn new(rng: &mut dyn rand::RngCore, world: HittableList, t0: f64, t1: f64) -> BVHNode {
        Self::new_impl(rng, world.objects.to_vec(), t0, t1)
    }

    fn new_impl(
        rng: &mut dyn rand::RngCore,
        mut list: Vec<Arc<dyn Hittable + Send + Sync>>,
        t0: f64,
        t1: f64,
    ) -> BVHNode {
        let uni = rand::distributions::Uniform::from(0..3);
        let axis = uni.sample(rng);
        let comparator = |a: &Arc<dyn Hittable + Send + Sync>,
                          b: &Arc<dyn Hittable + Send + Sync>| {
            let box_a = a.bounding_box(t0, t1).unwrap();
            let box_b = b.bounding_box(t0, t1).unwrap();
            match axis {
                0 => box_a.min.x.partial_cmp(&box_b.min.x).unwrap(),
                1 => box_a.min.y.partial_cmp(&box_b.min.y).unwrap(),
                2 => box_a.min.z.partial_cmp(&box_b.min.z).unwrap(),
                _ => panic!("invalid axis"),
            }
        };

        let (left, right) = if list.len() == 1 {
            (list[0].clone(), list[0].clone())
        } else if list.len() == 2 {
            if comparator(&list[0], &list[1]) == Ordering::Less {
                (list[0].clone(), list[1].clone())
            } else {
                (list[1].clone(), list[0].clone())
            }
        } else {
            list.sort_by(comparator);
            let mid = list.len() / 2;
            (
                Arc::new(BVHNode::new_impl(rng, list[..mid].to_vec(), t0, t1))
                    as Arc<dyn Hittable + Send + Sync>,
                Arc::new(BVHNode::new_impl(rng, list[mid..].to_vec(), t0, t1))
                    as Arc<dyn Hittable + Send + Sync>,
            )
        };

        let box_left = left.bounding_box(t0, t1).unwrap();
        let box_right = right.bounding_box(t0, t1).unwrap();
        let box_ = surrounding_box(&box_left, &box_right);

        BVHNode { left, right, box_ }
    }
}

impl Hittable for BVHNode {
    fn hit(&self, r: &Ray, t_min: f64, mut t_max: f64) -> Option<HitRecord> {
        if !self.box_.hit(r, t_min, t_max) {
            None
        } else {
            let left_rec = self.left.hit(r, t_min, t_max);
            if let Some(rec) = &left_rec {
                t_max = f64::min(t_max, rec.t);
            }

            self.right.hit(r, t_min, t_max).or_else(|| left_rec)
        }
    }
    fn bounding_box(&self, _t0: f64, _t1: f64) -> Option<AABB> {
        Some(self.box_.clone())
    }
}

#[derive(Default)]
pub struct HittableList {
    objects: Arc<Vec<Arc<dyn Hittable + Send + Sync>>>,
}

impl HittableList {
    pub fn new(objects: Vec<Arc<dyn Hittable + Send + Sync>>) -> HittableList {
        HittableList {
            objects: Arc::new(objects),
        }
    }
}

impl Hittable for HittableList {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut closest_so_far = t_max;
        let mut hit_record = None;

        for object in self.objects.iter() {
            if let Some(hit) = object.hit(r, t_min, closest_so_far) {
                closest_so_far = hit.t;
                hit_record = Some(hit);
            }
        }

        hit_record
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        self.objects
            .iter()
            .filter_map(|object| object.bounding_box(t0, t1))
            .reduce(|a, b| surrounding_box(&a, &b))
    }
}

pub struct XYRect {
    pub x0: f64,
    pub x1: f64,
    pub y0: f64,
    pub y1: f64,
    pub k: f64,
    pub material: Arc<dyn Material + Send + Sync>,
}

impl XYRect {
    pub fn new(
        x0: f64,
        x1: f64,
        y0: f64,
        y1: f64,
        k: f64,
        material: Arc<dyn Material + Send + Sync>,
    ) -> XYRect {
        XYRect {
            x0,
            x1,
            y0,
            y1,
            k,
            material,
        }
    }
}

impl Hittable for XYRect {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - r.orig.z) / r.dir.z;
        if t < t_min || t > t_max {
            return None;
        }

        let x = r.orig.x + t * r.dir.x;
        let y = r.orig.y + t * r.dir.y;

        if x < self.x0 || x > self.x1 || y < self.y0 || y > self.y1 {
            return None;
        }

        let u = (x - self.x0) / (self.x1 - self.x0);
        let v = (y - self.y0) / (self.y1 - self.y0);

        let outward_normal = Vec3::new(0.0, 0.0, 1.0);
        let p = r.at(t);

        Some(HitRecord::new(
            r,
            p,
            outward_normal,
            self.material.clone(),
            t,
            u,
            v,
        ))
    }

    fn bounding_box(&self, _t0: f64, _t1: f64) -> Option<AABB> {
        Some(AABB::new(
            Vec3::new(self.x0, self.y0, self.k - 0.0001),
            Vec3::new(self.x1, self.y1, self.k + 0.0001),
        ))
    }
}

pub struct YZRect {
    pub y0: f64,
    pub y1: f64,
    pub z0: f64,
    pub z1: f64,
    pub k: f64,
    pub material: Arc<dyn Material + Send + Sync>,
}

impl YZRect {
    pub fn new(
        y0: f64,
        y1: f64,
        z0: f64,
        z1: f64,
        k: f64,
        material: Arc<dyn Material + Send + Sync>,
    ) -> YZRect {
        YZRect {
            y0,
            y1,
            z0,
            z1,
            k,
            material,
        }
    }
}

impl Hittable for YZRect {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - r.orig.x) / r.dir.x;
        if t < t_min || t > t_max {
            return None;
        }

        let y = r.orig.y + t * r.dir.y;
        let z = r.orig.z + t * r.dir.z;

        if y < self.y0 || y > self.y1 || z < self.z0 || z > self.z1 {
            return None;
        }

        let u = (y - self.y0) / (self.y1 - self.y0);
        let v = (z - self.z0) / (self.z1 - self.z0);

        let outward_normal = Vec3::new(1.0, 0.0, 0.0);
        let p = r.at(t);

        Some(HitRecord::new(
            r,
            p,
            outward_normal,
            self.material.clone(),
            t,
            u,
            v,
        ))
    }

    fn bounding_box(&self, _t0: f64, _t1: f64) -> Option<AABB> {
        Some(AABB::new(
            Vec3::new(self.k - 0.0001, self.y0, self.z0),
            Vec3::new(self.k + 0.0001, self.y1, self.z1),
        ))
    }
}

pub struct ZXRect {
    pub z0: f64,
    pub z1: f64,
    pub x0: f64,
    pub x1: f64,
    pub k: f64,
    pub material: Arc<dyn Material + Send + Sync>,
}

impl ZXRect {
    pub fn new(
        z0: f64,
        z1: f64,
        x0: f64,
        x1: f64,
        k: f64,
        material: Arc<dyn Material + Send + Sync>,
    ) -> ZXRect {
        ZXRect {
            z0,
            z1,
            x0,
            x1,
            k,
            material,
        }
    }
}

impl Hittable for ZXRect {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - r.orig.y) / r.dir.y;
        if t < t_min || t > t_max {
            return None;
        }

        let z = r.orig.z + t * r.dir.z;
        let x = r.orig.x + t * r.dir.x;

        if z < self.z0 || z > self.z1 || x < self.x0 || x > self.x1 {
            return None;
        }

        let u = (z - self.z0) / (self.z1 - self.z0);
        let v = (x - self.x0) / (self.x1 - self.x0);

        let outward_normal = Vec3::new(0.0, 1.0, 0.0);
        let p = r.at(t);

        Some(HitRecord::new(
            r,
            p,
            outward_normal,
            self.material.clone(),
            t,
            u,
            v,
        ))
    }

    fn bounding_box(&self, _t0: f64, _t1: f64) -> Option<AABB> {
        Some(AABB::new(
            Vec3::new(self.x0, self.k - 0.0001, self.z0),
            Vec3::new(self.x1, self.k + 0.0001, self.z1),
        ))
    }
}

pub struct Box {
    pub box_min: Vec3,
    pub box_max: Vec3,
    pub sides: HittableList,
}

impl Box {
    pub fn new(box_min: Vec3, box_max: Vec3, material: Arc<dyn Material + Send + Sync>) -> Box {
        let mut sides: Vec<Arc<dyn Hittable + Send + Sync>> = Vec::new();
        sides.push(Arc::new(XYRect::new(
            box_min.x,
            box_max.x,
            box_min.y,
            box_max.y,
            box_min.z,
            material.clone(),
        )));
        sides.push(Arc::new(XYRect::new(
            box_min.x,
            box_max.x,
            box_min.y,
            box_max.y,
            box_max.z,
            material.clone(),
        )));
        sides.push(Arc::new(YZRect::new(
            box_min.y,
            box_max.y,
            box_min.z,
            box_max.z,
            box_min.x,
            material.clone(),
        )));
        sides.push(Arc::new(YZRect::new(
            box_min.y,
            box_max.y,
            box_min.z,
            box_max.z,
            box_max.x,
            material.clone(),
        )));
        sides.push(Arc::new(ZXRect::new(
            box_min.z,
            box_max.z,
            box_min.x,
            box_max.x,
            box_min.y,
            material.clone(),
        )));
        sides.push(Arc::new(ZXRect::new(
            box_min.z,
            box_max.z,
            box_min.x,
            box_max.x,
            box_max.y,
            material.clone(),
        )));

        Box {
            box_min,
            box_max,
            sides: HittableList::new(sides),
        }
    }
}

impl Hittable for Box {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        self.sides.hit(r, t_min, t_max)
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        Some(AABB::new(self.box_min, self.box_max))
    }
}

pub struct Translate {
    pub obj: Arc<dyn Hittable + Send + Sync>,
    pub offset: Vec3,
}

impl Translate {
    pub fn new(obj: Arc<dyn Hittable + Send + Sync>, offset: Vec3) -> Translate {
        Translate { obj, offset }
    }
}

impl Hittable for Translate {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let moved_r = Ray::new(r.orig - self.offset, r.dir, r.time);
        let rec = self.obj.hit(&moved_r, t_min, t_max)?;

        Some(HitRecord::new(
            &moved_r,
            rec.p + self.offset,
            rec.normal,
            rec.material,
            rec.t,
            rec.u,
            rec.v,
        ))
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        self.obj
            .bounding_box(t0, t1)
            .map(|aabb| AABB::new(aabb.min + self.offset, aabb.max + self.offset))
    }
}

pub struct RotateY {
    pub ptr: Arc<dyn Hittable + Send + Sync>,
    pub sin_theta: f64,
    pub cos_theta: f64,
}

impl RotateY {
    pub fn new(ptr: Arc<dyn Hittable + Send + Sync>, angle: f64) -> RotateY {
        let radians = angle.to_radians();
        RotateY {
            ptr,
            sin_theta: radians.sin(),
            cos_theta: radians.cos(),
        }
    }
}

impl Hittable for RotateY {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let origin = Vec3::new(
            self.cos_theta * r.orig.x - self.sin_theta * r.orig.z,
            r.orig.y,
            self.sin_theta * r.orig.x + self.cos_theta * r.orig.z,
        );

        let direction = Vec3::new(
            self.cos_theta * r.dir.x - self.sin_theta * r.dir.z,
            r.dir.y,
            self.sin_theta * r.dir.x + self.cos_theta * r.dir.z,
        );
        let r = Ray::new(origin, direction, r.time);
        let rec = self.ptr.hit(&r, t_min, t_max)?;

        let p = Vec3::new(
            self.cos_theta * rec.p.x + self.sin_theta * rec.p.z,
            rec.p.y,
            -self.sin_theta * rec.p.x + self.cos_theta * rec.p.z,
        );

        let normal = Vec3::new(
            self.cos_theta * rec.normal.x + self.sin_theta * rec.normal.z,
            rec.normal.y,
            -self.sin_theta * rec.normal.x + self.cos_theta * rec.normal.z,
        );

        Some(HitRecord::new(
            &r,
            p,
            normal,
            rec.material,
            rec.t,
            rec.u,
            rec.v,
        ))
    }

    fn bounding_box(&self, t0: f64, t1: f64) -> Option<AABB> {
        let bbox = self.ptr.bounding_box(t0, t1)?;

        let v = iproduct!(
            [bbox.min.x, bbox.max.x],
            [bbox.min.y, bbox.max.y],
            [bbox.min.z, bbox.max.z]
        )
        .map(|(x, y, z)| {
            Vec3::new(
                self.cos_theta * x + self.sin_theta * z,
                y,
                -self.sin_theta * x + self.cos_theta * z,
            )
        })
        .collect::<Vec<_>>();
        let min = Vec3::new(
            v.iter().map(|v| v.x).reduce(|a, b| a.min(b))?,
            v.iter().map(|v| v.y).reduce(|a, b| a.min(b))?,
            v.iter().map(|v| v.z).reduce(|a, b| a.min(b))?,
        );
        let max = Vec3::new(
            v.iter().map(|v| v.x).reduce(|a, b| a.max(b))?,
            v.iter().map(|v| v.y).reduce(|a, b| a.max(b))?,
            v.iter().map(|v| v.z).reduce(|a, b| a.max(b))?,
        );

        Some(AABB::new(min, max))
    }
}
