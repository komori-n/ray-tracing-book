const WIDTH: usize = 256;
const HEIGHT: usize = 256;

mod raytrace;

use raytrace::color::{ppm_string, Color};

fn output() {
    println!("P3");
    println!("{} {}", WIDTH, HEIGHT);
    println!("255");

    for j in (0..HEIGHT).rev() {
        eprint!("\rScanlines remaining: {}   ", j);
        for i in 0..WIDTH {
            let r = (i as f64) / ((WIDTH - 1) as f64);
            let g = (j as f64) / ((HEIGHT - 1) as f64);
            let b = 0.25;
            let c = Color::new(r, g, b);

            println!("{}", ppm_string(c));
        }
    }
    eprintln!("");
    eprintln!("Done");
}

fn main() {
    output();
}
