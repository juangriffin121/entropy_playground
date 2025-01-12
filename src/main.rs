use bases::{plot::plot_recipient, recipient::Recipient};
use rand::{rngs::StdRng, SeedableRng};

#[allow(unused)]
mod bases;
fn main() {
    //something fishy about rngs
    let seed = 42; // Replace with a desired seed value for reproducibility
    let mut rng1 = StdRng::seed_from_u64(seed);
    let mut recipient = Recipient::from_json("simulation_data/recipient.json", &mut rng1);

    for i in 0..recipient.iterations {
        println!("bacwards {:?} {:?}", i, recipient.contents.len());
        recipient.update(recipient.dt);
        if i % recipient.simulation_rate == 0 {
            plot_recipient(
                &recipient,
                &format!(
                    "plots/{:05}.png",
                    recipient.iterations / recipient.simulation_rate
                        - i / recipient.simulation_rate
                ),
                true,
            );
        }
    }

    let mut rng2 = StdRng::seed_from_u64(seed);
    let mut rec2 = Recipient::from_json("simulation_data/recipient.json", &mut rng2);
    rec2.dt = -rec2.dt;
    for i in 1..rec2.iterations {
        println!("forwards {:?} {:?}", i, rec2.contents.len());
        rec2.update(rec2.dt);
        if i % rec2.simulation_rate == 0 {
            plot_recipient(
                &rec2,
                &format!(
                    "plots/{:05}.png",
                    recipient.iterations / rec2.simulation_rate + i / rec2.simulation_rate
                ),
                true,
            );
        }
    }
}
