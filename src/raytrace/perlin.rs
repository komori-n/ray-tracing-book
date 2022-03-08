use rand::prelude::{Distribution, SliceRandom};

use super::vec3::{dot, Vec3};

const POINT_COUNT: usize = 256;

pub struct Perlin {
    ranvec: Vec<Vec3>,
    perm_x: Vec<usize>,
    perm_y: Vec<usize>,
    perm_z: Vec<usize>,
}

impl Perlin {
    pub fn new(rng: &mut dyn rand::RngCore) -> Perlin {
        let uni = rand::distributions::Uniform::from(0.0..1.0);
        let ranvec = (0..POINT_COUNT)
            .map(|_| Vec3::random_with_minmax(rng, -1.0, 1.0))
            .collect::<Vec<_>>();
        let perm_x = Perlin::generate_perm(rng);
        let perm_y = Perlin::generate_perm(rng);
        let perm_z = Perlin::generate_perm(rng);

        Perlin {
            ranvec,
            perm_x,
            perm_y,
            perm_z,
        }
    }

    pub fn noise(&self, p: &Vec3) -> f64 {
        let u = p.x - p.x.floor();
        let v = p.y - p.y.floor();
        let w = p.z - p.z.floor();

        let uu = u * u * (3.0 - 2.0 * u);
        let vv = v * v * (3.0 - 2.0 * v);
        let ww = w * w * (3.0 - 2.0 * w);

        let i = (p.x.floor() as i64).to_le_bytes()[0] as usize;
        let j = (p.y.floor() as i64).to_le_bytes()[0] as usize;
        let k = (p.z.floor() as i64).to_le_bytes()[0] as usize;

        let mut accum = 0.0;
        for ii in 0..2 {
            for jj in 0..2 {
                for kk in 0..2 {
                    let ti = i + ii;
                    let tj = j + jj;
                    let tk = k + kk;
                    let idx = self.perm_x[ti % 256] ^ self.perm_y[tj % 256] ^ self.perm_z[tk % 256];
                    let val = self.ranvec[idx];
                    accum += ((ii as f64) * uu + (1.0 - (ii as f64)) * (1.0 - uu))
                        * ((jj as f64) * vv + (1.0 - (jj as f64)) * (1.0 - vv))
                        * ((kk as f64) * ww + (1.0 - (kk as f64)) * (1.0 - ww))
                        * dot(
                            val,
                            Vec3::new(u - (ii as f64), v - (jj as f64), w - (kk as f64)),
                        );
                }
            }
        }

        accum
    }

    pub fn turb(&self, p: &Vec3, depth: i32) -> f64 {
        let (acc, _, _) = (0..depth).fold((0.0, 1.0, *p), |(acc, weight, tmp_p), _| {
            (acc + weight * self.noise(&tmp_p), weight * 0.5, tmp_p * 2.0)
        });

        acc.abs()
    }

    fn generate_perm(rng: &mut dyn rand::RngCore) -> Vec<usize> {
        let mut v = (0..POINT_COUNT).collect::<Vec<_>>();
        v.shuffle(rng);
        v
    }
}
