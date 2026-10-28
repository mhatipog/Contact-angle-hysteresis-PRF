"""Optional research figure from the reproducible reduced volume-cycle model.

pip install numpy scipy matplotlib
python research/microscopic_wetting/plot_cycle.py --output-prefix results/cah
"""
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from quasistatic_cycle import volume_cycle


def make_plots(prefix, steps=544):
    rows,jumps,work,dE,diss=volume_cycle(steps_per_leg=steps)
    midpoint=steps+1
    q=(rows[:,1]+rows[:,2])/2
    p=rows[:,4]
    prefix=Path(prefix)
    prefix.parent.mkdir(parents=True,exist_ok=True)

    fig,ax=plt.subplots(figsize=(8.4,5.4))
    ax.plot(rows[:midpoint,0],q[:midpoint],label="Volume decreasing")
    ax.plot(rows[midpoint-1:,0],q[midpoint-1:],label="Volume increasing")
    for x in (.213851157109733,.2076472661283195):
        ax.axvline(x,linestyle=":",alpha=.5)
    ax.scatter([rows[0,0],rows[-1,0]],[q[0],q[-1]],zorder=5)
    ax.set(xlabel="Dimensionless liquid volume V",
           ylabel="Droplet center (a+b)/2",
           title="Two sequential release/repinning events and one-way relocation")
    ax.grid(alpha=.25);ax.legend()
    fig.tight_layout()
    fig.savefig(str(prefix)+"_center.png",dpi=180)
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.4,5.4))
    ax.plot(rows[:midpoint,0],p[:midpoint],label="Volume decreasing")
    ax.plot(rows[midpoint-1:,0],p[midpoint-1:],label="Volume increasing")
    ax.set(xlabel="Dimensionless liquid volume V",
           ylabel="Conjugate pressure p=dE/dV",
           title="Metastable pressure-volume path")
    ax.grid(alpha=.25);ax.legend()
    fig.tight_layout()
    fig.savefig(str(prefix)+"_pressure.png",dpi=180)
    plt.close(fig)

    print(f"{len(jumps)} repinning events; irreversible energy loss {diss:.12g}")
    return jumps


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-prefix",default="cah")
    parser.add_argument("--steps",type=int,default=544)
    args=parser.parse_args()
    make_plots(args.output_prefix,args.steps)
