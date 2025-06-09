# Hot-Fibrosis Fixed-Point Scan  
_Reproduction & extension of Adler et al. (2020) wound-healing circuit_

This repository contains **one self-contained script**  
`lambda_scan_standalone.py`  
that sweeps the two key feedback parameters of the Adler-2020 model

* **λ₁** – PDGF autocrine strength (fibroblast self-stimulation)  
* **λ₂** – CSF-1 mediated macrophage proliferation

and isolates the **hot-fibrosis** fixed point (highest fibroblast burden **F\***),
reporting the corresponding macrophage load **M\*** and exporting a publication-ready
heat-map:

![example output](hot_lambda_dependence.png)

---

## Quick start

```bash
git clone git@github.com:Tom-Brenner/hot-fibrosis-scan.git
cd hot-fibrosis-scan
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt      # numpy scipy matplotlib
python lambda_scan_standalone.py \
       --lambda1_lo 0.2  --lambda1_hi 2 \
       --lambda2_lo 0.2  --lambda2_hi 2 \
       --N 30 --nproc 8   > scan.csv
