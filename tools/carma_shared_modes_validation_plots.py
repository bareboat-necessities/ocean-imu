"""Compact vector figures and full-rate-score tables for the CARMA article."""
from __future__ import annotations

import json
import hashlib
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "plots"))
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family": "serif", "font.size": 8,
    "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "axes.grid": True, "grid.alpha": .23,
    "axes.unicode_minus": False, "pgf.texsystem": "pdflatex", "pgf.rcfonts": False})
import matplotlib.pyplot as plt
from svg_determinism import configure_svg
from vector_export import save_figure
configure_svg(matplotlib)

LABEL = {"CSM": "Shared CARMA", "Independent": "Uncoupled CARMA", "OU-II": "OU-II",
         "OU-III": "OU-III", "TFG": "TFG", "PII": "PII", "TVG-NLO": "TVG-NLO"}
TEX = {**LABEL, "OU-II": "OU--II", "OU-III": "OU--III", "TVG-NLO": "TVG--NLO"}
COLOR = dict(zip(LABEL, ("#0072b2", "#777777", "#009e73", "#d55e00", "#cc79a7", "#e69f00", "#56b4e9")))
KALMAN = ("CSM", "Independent", "OU-II", "OU-III", "TFG")
METRICS = ("disp_z_pct_hs", "disp_3d_pct_hs", "pitch_rms_deg")


def trace(results, name):
    return pd.read_csv(results / "traces" / (name + ".csv.gz"))


def save(fig, out, name):
    fig.align_ylabels()
    fig.tight_layout(pad=.65)
    for ext in ("svg", "pgf"):
        kwargs = {"metadata": {"Date": None}} if ext == "svg" else {}
        target=out/(name+"."+ext)
        temporary=out/(name+"."+ext+".pending")
        # Readers (including an article build) must never see a partial PGF.
        with temporary.open("w",encoding="utf-8") as stream:
            save_figure(fig, stream, format=ext, bbox_inches="tight", **kwargs)
        ending=temporary.read_text().rstrip()
        if not ending.endswith(r"\endgroup%" if ext=="pgf" else "</svg>"):
            raise RuntimeError(f"incomplete {ext} export: {target}")
        if ext == "svg":
            temporary.write_text("\n".join(line.rstrip() for line in ending.splitlines())+"\n")
        temporary.replace(target)
    plt.close(fig)


def paired(ax, data, columns, ylabel):
    for n, (ref, est) in enumerate(columns):
        color = "black" if len(columns) == 1 else ("#0072b2", "#d55e00")[n]
        suffix = "" if len(columns) == 1 else (" X", " Y")[n]
        ax.plot(data.time, data[ref], color=color, lw=.9, label="Truth" + suffix)
        ax.plot(data.time, data[est], color=(COLOR["CSM"] if len(columns) == 1 else color),
                ls="--", lw=.8, label="Estimate" + suffix)
    ax.set_ylabel(ylabel)
    ax.legend(loc="upper right", ncol=2)


def nominal_figures(results, out):
    for axes_name in ("xy", "z"):
        fig, axes = plt.subplots(3,1,figsize=(3.45,4.05),sharex=True)
        for ax, kind, unit in zip(axes,("disp","vel","acc"),("m","m/s",r"m/s$^2$")):
            d=trace(results,f"nominal_CSM_{kind}_{axes_name}")
            pairs=[(f"{kind}_ref_{a}",f"{kind}_est_{a}") for a in axes_name]
            paired(ax,d,pairs,dict(disp="Position",vel="Velocity",acc="Acceleration")[kind]+f" [{unit}]")
        axes[-1].set_xlabel("Time [s]")
        save(fig,out,"csm_validation_"+axes_name)
    fig, axes=plt.subplots(3,1,figsize=(3.45,3.60),sharex=True)
    for ax, axis in zip(axes,("roll","pitch","yaw")):
        paired(ax,trace(results,"nominal_CSM_"+axis),[(axis+"_ref",axis+"_est")],axis.capitalize()+" [deg]")
    axes[-1].set_xlabel("Time [s]")
    save(fig,out,"csm_validation_attitude")
    for kind,unit in (("acc",r"m/s$^2$"),("gyro","rad/s")):
        fig, axes=plt.subplots(3,1,figsize=(3.45,3.60),sharex=True)
        for ax,axis in zip(axes,"xyz"):
            d=trace(results,f"nominal_CSM_{kind}_bias_{axis}")
            paired(ax,d,[(f"{kind}_bias_{axis}",f"{kind}_bias_est_{axis}")],axis.upper()+f" bias [{unit}]")
        axes[-1].set_xlabel("Time [s]")
        save(fig,out,"csm_validation_"+kind+"_bias")
    fig,axes=plt.subplots(3,1,figsize=(3.45,3.60),sharex=True)
    for ax,col,label in zip(axes,("Tp_tuner_s","sigma_a_applied","R_p0_applied"),
             ("Carrier period [s]",r"Applied $\sigma_a$ [m/s$^2$]",r"Frequency scale $\kappa$")):
        d=trace(results,"nominal_CSM_"+col)
        ax.plot(d.time,d[col],color=COLOR["CSM"],lw=.9)
        ax.set_ylabel(label)
    axes[-1].set_xlabel("Time [s]")
    save(fig,out,"csm_validation_tuning")
    fig,axes=plt.subplots(2,1,figsize=(6.9,3.4),sharex=True)
    groups=(("CSM","Independent","OU-III"),("OU-II","TFG","PII","TVG-NLO"))
    for ax,group in zip(axes,groups):
        ref=trace(results,"nominal_CSM_z")
        ax.plot(ref.time,ref.disp_ref_z,color="black",lw=1.2,label="Truth")
        for method in group:
            d=trace(results,f"nominal_{method}_z")
            ax.plot(d.time,d.disp_est_z,color=COLOR[method],lw=.8,ls="--" if method=="CSM" else "-",label=LABEL[method])
        ax.set_ylabel("Vertical position [m]")
        ax.legend(ncol=3,loc="upper right")
    axes[-1].set_xlabel("Time [s]")
    save(fig,out,"csm_validation_observers")


def mean_interval(values):
    a=np.asarray(values,dtype=float)
    if len(a)<2 or not np.isfinite(a).all():
        raise ValueError("paired interval requires complete finite seed-level data")
    rng=np.random.default_rng(20261108)
    boot=a[rng.integers(0,len(a),(10000,len(a)))].mean(axis=1)
    return float(a.mean()),float(np.percentile(boot,2.5)),float(np.percentile(boot,97.5))


def comparative_figure(df,out):
    d=df[(df.study=="stationary") & df.method.isin(KALMAN)]
    fig,axes=plt.subplots(3,2,figsize=(6.9,5.55),sharex=True)
    for j,spectrum in enumerate(("JONSWAP","PM-Stokes")):
        for i,(metric,ylabel) in enumerate(zip(METRICS,(r"Heave error / $H_s$ [%]",r"3-D error / $H_s$ [%]","Pitch RMS [deg]"))):
            ax=axes[i,j]
            for method in KALMAN:
                g=d[(d.spectrum==spectrum)&(d.method==method)].groupby("hs_m")[metric]
                means=g.mean()
                intervals=np.array([mean_interval(g.get_group(h))[1:] for h in means.index])
                ax.errorbar(means.index,means,yerr=np.vstack((means.values-intervals[:,0],intervals[:,1]-means.values)),
                            color=COLOR[method],lw=.9,marker="o",ms=2.5,capsize=2,label=LABEL[method])
            ax.set_xscale("log")
            ax.set_xticks([.27,1.5,4,8.5],["0.27","1.5","4","8.5"])
            if j==0: ax.set_ylabel(ylabel)
            if i==0: ax.set_title(spectrum,fontsize=9)
            if i==2: ax.set_xlabel(r"Incident significant height $H_s$ [m]")
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.5,1.045),ncol=3,frameon=False)
    save(fig,out,"csm_validation_comparison")


def pooled(values):
    a=np.asarray(values,dtype=float)
    return float(np.sqrt(np.mean(a*a)))


def engine_summary(df):
    rows=[]
    for (method,guard,condition),g in df[df.study=="engine"].groupby(["method","guard","condition"]):
        if len(g)!=8:
            raise ValueError(f"incomplete engine cell {method}/{guard}/{condition}")
        row={"method":method,"guard":guard,"condition":condition,"records":len(g)}
        for name in ("disp_z_rms_m","disp_3d_rms_m","roll_rms_deg","pitch_rms_deg","yaw_rms_deg"):
            row[name]=pooled(g[name])
        for name in ("rpm","level_mps2","bandwidth_hz","recorded_rms_mps2"):
            row[name]=float(g[name].iloc[0])
        row["disp_offset_m"]=float(np.sqrt(np.mean(sum(g[f"disp_{a}_mean_m"]**2 for a in "xyz"))))
        row["pitch_offset_deg"]=pooled(g.pitch_mean_deg)
        rows.append(row)
    return pd.DataFrame(rows)


def engine_figures(e,out):
    fig,axes=plt.subplots(3,1,figsize=(3.45,4.6),sharex=True)
    for method,guard,style in (("CSM","off","-"),("Independent","off","-"),("OU-III","off","-"),("CSM","guard","--"),("CSM","deployed",":")):
        d=e[(e.method==method)&(e.guard==guard)&e.condition.isin(["rpm800","rpm1600","rpm2400","rpm3200"])].sort_values("rpm")
        label=LABEL[method]+({"off":" (raw)","guard":" (guard)","deployed":" (guard + R)"}[guard])
        for ax,metric,ylabel in zip(axes,("disp_3d_rms_m","pitch_rms_deg","yaw_rms_deg"),("3-D RMS [m]","Pitch RMS [deg]","Yaw RMS [deg]")):
            ax.plot(d.rpm,d[metric],color=COLOR[method],ls=style,marker="o",ms=2.5,lw=1,label=label)
            base=e[(e.method==method)&(e.guard==guard)&(e.condition=="off")]
            if guard=="off": ax.axhline(base[metric].iloc[0],color=COLOR[method],ls=":",lw=.6,alpha=.65)
            ax.set_ylabel(ylabel)
            ax.set_yscale("log")
    axes[0].legend(loc="upper left",fontsize=6)
    axes[-1].set_xlabel("Engine speed [rpm]")
    save(fig,out,"csm_validation_engine_speed")
    fig,axes=plt.subplots(3,1,figsize=(3.45,4.6))
    panels=(("level0.3","rpm2400","level1.2"),("bw20","bw40","rpm2400","bw160"),("matched20","matched40","rpm2400","matched160"))
    for ax,conditions,xname,title in zip(axes,panels,("level_mps2","bandwidth_hz","bandwidth_hz"),("Hull level","Sensor bandwidth","Matched recorded power")):
        for method in ("CSM","Independent","OU-III"):
            d=e[(e.method==method)&(e.guard=="off")&e.condition.isin(conditions)].sort_values(xname)
            ax.plot(d[xname],d.disp_3d_rms_m,color=COLOR[method],marker="o",ms=2.5,label=LABEL[method],lw=.9)
        ax.set_yscale("log")
        ax.set_title(title,fontsize=8)
        ax.set_ylabel("3-D RMS [m]")
        ax.set_xlabel(r"Hull RMS [m/s$^2$]" if xname=="level_mps2" else "Bandwidth [Hz]")
    axes[0].set_ylabel("3-D RMS [m]")
    axes[0].legend(fontsize=6)
    save(fig,out,"csm_validation_engine_mechanism")


def transition_figure(results,out):
    import ou_roundtrip_transition as rt
    import ou_validation as ov
    fig,axes=plt.subplots(3,2,figsize=(6.9,5.85),sharex=True)
    t=np.linspace(0,1200,2401)
    weight,_,_=rt.roundtrip_profile(ov,t)
    height=np.array([ov.mixture_significant_height_m(1.5,4.0,float(w)) for w in weight])
    axes[0,0].plot(t,height,color="black",lw=1)
    axes[0,0].set_ylabel(r"Mixture $H_s$ [m]")
    d=trace(results,"roundtrip_CSM_z")
    paired(axes[0,1],d,[("disp_ref_z","disp_est_z")],"Heave [m]")
    for method in ("CSM","Independent","OU-III"):
        d=trace(results,f"roundtrip_{method}_z")
        axes[1,0].plot(d.time,d.disp_est_z-d.disp_ref_z,color=COLOR[method],lw=.6,label=LABEL[method])
    axes[1,0].set_ylabel("Heave error [m]")
    axes[1,0].legend(fontsize=6)
    d=trace(results,"roundtrip_CSM_R_p0_applied")
    axes[1,1].plot(d.time,d.R_p0_applied,color=COLOR["CSM"],lw=.9)
    axes[1,1].set_ylabel(r"Applied $\kappa$")
    for col,color in (("roll","#0072b2"),("pitch","#d55e00"),("yaw","#009e73")):
        d=trace(results,"roundtrip_CSM_"+col)
        error=(d[col+"_est"]-d[col+"_ref"]+180)%360-180
        axes[2,0].plot(d.time,error,color=color,lw=.65,label=col)
    axes[2,0].set_ylabel("Attitude error [deg]")
    axes[2,0].legend(fontsize=6,ncol=3)
    d=trace(results,"roundtrip_CSM_sigma_a_applied")
    axes[2,1].plot(d.time,d.sigma_a_applied,color=COLOR["CSM"],lw=.9)
    axes[2,1].set_ylabel(r"Applied $\sigma_a$ [m/s$^2$]")
    for ax in axes.flat:
        for a,b in ((400,520),(800,920)):
            ax.axvspan(a,b,color="0.8",alpha=.4,zorder=0)
        ax.axvline(300,color="0.5",ls="-.",lw=.6)
        ax.set_xlim(0,1200)
    for ax in axes[-1]: ax.set_xlabel("Time [s]")
    save(fig,out,"csm_validation_roundtrip")


def pm(values):
    a=np.asarray(values,dtype=float)
    return f"${a.mean():.3f}\\pm{a.std(ddof=1):.3f}$"


def tables(df,e,out):
    lines=["% Generated by tools/carma_shared_modes_validation.py; full-rate scores."]
    s=df[df.study=="stationary"]
    nominal=s[(s.method=="CSM")&(s.spectrum=="JONSWAP")&(s.hs_m==1.5)]
    lines += [r"\begin{table}[t]\centering\footnotesize\setlength{\tabcolsep}{3pt}",
        r"\caption{Shared CARMA on medium JONSWAP: ten-seed mean $\pm$ sample standard deviation of RMS over the final 900 s.}",
        r"\label{tab:csm-validation-nominal}",r"\begin{tabular}{@{}llr@{}}\toprule",
        r"Quantity & Unit & RMS error\\\midrule"]
    for key,label,unit in (("disp_x_rms_m","X displacement","m"),("disp_y_rms_m","Y displacement","m"),("disp_z_rms_m","Heave","m"),
            ("roll_rms_deg","Roll","deg"),("pitch_rms_deg","Pitch","deg"),("yaw_rms_deg","Yaw","deg"),
            ("accel_bias_3d_rms_mps2","Accel. bias (3-D)",r"m/s$^2$"),("gyro_bias_3d_rms_radps","Gyro bias (3-D)","rad/s")):
        value=pm(nominal[key]) if key!="gyro_bias_3d_rms_radps" else f"$({nominal[key].mean()*1e4:.2f}\\pm{nominal[key].std(ddof=1)*1e4:.2f})\\,10^{{-4}}$"
        lines.append(f"{label} & {unit} & {value}"+r"\\")
    lines += [r"\bottomrule\end{tabular}\end{table}"]
    lines += [r"\begin{table*}[t]\centering\footnotesize",
        r"\caption{Paired ten-seed stationary replays. Mean $\pm$ sample standard deviation of the seed-level average over the four heights in each spectral family. Displacement errors are percentages of incident $H_s$; pitch is in degrees.}",
        r"\label{tab:csm-validation-comparison}",r"\begin{tabular}{@{}lrrrrrr@{}}\toprule",
        r"&\multicolumn{3}{c}{JONSWAP}&\multicolumn{3}{c}{PM--Stokes}\\",
        r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",r"Method & Heave [\%] & 3-D [\%] & Pitch [deg] & Heave [\%] & 3-D [\%] & Pitch [deg]\\\midrule"]
    for method in KALMAN:
        values=[]
        for spectrum in ("JONSWAP","PM-Stokes"):
            g=s[(s.method==method)&(s.spectrum==spectrum)].groupby("wave_seed")[list(METRICS)].mean()
            values += [pm(g[m]) for m in METRICS]
        lines.append(TEX[method]+" & "+" & ".join(values)+r"\\")
    lines += [r"\bottomrule\end{tabular}\end{table*}"]
    effects=[]
    lines += [r"\begin{table}[t]\centering\footnotesize\setlength{\tabcolsep}{2pt}",
        r"\caption{Shared CARMA minus comparator: seed-level mean differences and paired percentile-bootstrap 95\% intervals. Units are percentage points of $H_s$; uncoupled denotes the uncoupled CARMA arm.}",
        r"\label{tab:csm-validation-paired}",r"\begin{tabular}{@{}lllr@{}}\toprule",
        r"Spectrum & Comparator & Error & Difference [95\% CI]\\\midrule"]
    for spectrum in ("JONSWAP","PM-Stokes"):
        for comparator in ("OU-III","Independent"):
            for metric,label in (("disp_z_pct_hs","Heave"),("disp_3d_pct_hs","3-D")):
                g=s[(s.spectrum==spectrum)&s.method.isin(["CSM",comparator])].groupby(["wave_seed","method"])[metric].mean().unstack()
                mean,lo,hi=mean_interval(g.CSM-g[comparator])
                effects.append(dict(spectrum=spectrum,comparator=comparator,endpoint=metric,mean=mean,ci_low=lo,ci_high=hi))
                comp="OU--III" if comparator=="OU-III" else "Uncoupled"
                lines.append(f"{spectrum.replace('PM-Stokes','PM--Stokes')} & {comp} & {label} & ${mean:.3f}\\,[{lo:.3f},{hi:.3f}]$"+r"\\")
    lines += [r"\bottomrule\end{tabular}\end{table}"]
    lines += [r"\begin{table}[t]\centering\footnotesize\setlength{\tabcolsep}{2pt}",
        r"\caption{Common-output single-draw comparisons (wave/IMU/initialization seeds 11/101/1009). Arithmetic mean over four heights; heave as a percentage of incident $H_s$, pitch in degrees. These rows are not ten-seed estimates.}",
        r"\label{tab:csm-validation-observers}",r"\begin{tabular}{@{}lrrrr@{}}\toprule",
        r"&\multicolumn{2}{c}{JONSWAP}&\multicolumn{2}{c}{PM--Stokes}\\",
        r"Method & $Z$ [\%] & Pitch & $Z$ [\%] & Pitch\\\midrule"]
    for method in LABEL:
        values=[]
        for spectrum in ("JONSWAP","PM-Stokes"):
            g=s[(s.method==method)&(s.spectrum==spectrum)&(s.wave_seed==11)]
            values += [f"{g[m].mean():.3f}" for m in ("disp_z_pct_hs","pitch_rms_deg")]
        lines.append(TEX[method]+" & "+" & ".join(values)+r"\\")
    lines += [r"\bottomrule\end{tabular}\end{table}"]
    lines += [r"\begin{table*}[t]\centering\footnotesize",
        r"\caption{Nominal engine condition (2400 rpm, hull RMS 0.60 m/s$^2$, sensor bandwidth 80 Hz) against engine off. Pooled RMS over eight records, one common seed triplet. Offset is the RMS across records of the norm of each record's mean displacement error.}",
        r"\label{tab:csm-validation-engine}",r"\begin{tabular}{@{}llrrrrr@{}}\toprule",
        r"Method & Input / mitigation & Heave [m] & 3-D [m] & Offset [m] & Pitch [deg] & Yaw [deg]\\\midrule"]
    for method in ("CSM","Independent","OU-III"):
        for condition,guard,label in (("off","off","Engine off"),("rpm2400","off","Raw vibration"),
                ("accel-only","off","Accelerometer vibration only"),("rpm2400","guard","Vibration guard"),
                ("rpm2400","deployed","Guard + covariance")):
            d=e[(e.method==method)&(e.condition==condition)&(e.guard==guard)]
            if d.empty: continue
            r=d.iloc[0]
            values=[f"{r[m]:.3f}" for m in ("disp_z_rms_m","disp_3d_rms_m","disp_offset_m","pitch_rms_deg","yaw_rms_deg")]
            lines.append(TEX[method]+" & "+label+" & "+" & ".join(values)+r"\\")
    lines += [r"\bottomrule\end{tabular}\end{table*}"]
    lines += [r"\begin{table*}[t]\centering\footnotesize",
        r"\caption{Low--high--low transition: ten-seed mean $\pm$ sample standard deviation of full-rate RMS. Rise and fall, recovery and settled intervals remain separate. All errors are in metres.}",
        r"\label{tab:csm-validation-roundtrip}",r"\begin{tabular}{@{}lrrrrrr@{}}\toprule",
        r"&\multicolumn{2}{c}{Shared CARMA}&\multicolumn{2}{c}{Uncoupled CARMA}&\multicolumn{2}{c}{OU--III}\\",
        r"Interval [s] & Heave & 3-D & Heave & 3-D & Heave & 3-D\\\midrule"]
    import ou_roundtrip_transition as rt
    for segment,start,stop in rt.roundtrip_segments():
        values=[]
        for method in ("CSM","Independent","OU-III"):
            g=df[(df.study=="roundtrip")&(df.method==method)]
            values += [pm(g[f"seg_{segment}_{metric}"]) for metric in ("disp_z_rms_m","disp_3d_rms_m")]
        label=rt.SEGMENT_LABELS[segment]+f" ({start:g}--{stop:g})"
        lines.append(label+" & "+" & ".join(values)+r"\\")
    lines += [r"\bottomrule\end{tabular}\end{table*}"]
    blocks=re.findall(r"\\begin\{table\*?\}.*?\\end\{table\*?\}","\n".join(lines),flags=re.S)
    names=("Nominal","Comparison","Paired","Observers","Engine","Roundtrip")
    if len(blocks)!=len(names):
        raise RuntimeError("unexpected publication table set")
    contents=[lines[0]]
    for name,block in zip(names,blocks):
        contents.append("\\newcommand{\\CSMValidation"+name+"Table}{%\n"+block+"\n}")
    for effect in effects:
        spectrum="Jonswap" if effect["spectrum"]=="JONSWAP" else "PmStokes"
        endpoint="Heave" if effect["endpoint"]=="disp_z_pct_hs" else "ThreeD"
        name="CSM"+spectrum+endpoint+("Uncoupled" if effect["comparator"]=="Independent" else "")
        contents.append("\\newcommand{\\"+name+"Difference}{"+f"{effect['mean']:.3f}"+"}")
        contents.append("\\newcommand{\\"+name+"Interval}{"+f"[{effect['ci_low']:.3f},{effect['ci_high']:.3f}]"+"}")
    def engine_value(condition,guard,field):
        return float(e[(e.method=="CSM")&(e.condition==condition)&(e.guard==guard)][field].iloc[0])
    baseline=engine_value("off","off","disp_3d_rms_m")
    for name,value in (("CSMEngineBaselineThreeD",baseline),
          ("CSMEngineRawThreeD",engine_value("rpm2400","off","disp_3d_rms_m")),
          ("CSMEngineGuardOnlyThreeD",engine_value("rpm2400","guard","disp_3d_rms_m")),
          ("CSMEngineDeployedThreeD",engine_value("rpm2400","deployed","disp_3d_rms_m")),
          ("CSMEngineGyroPathChange",100*abs(engine_value("accel-only","off","disp_3d_rms_m")/
                                          engine_value("rpm2400","off","disp_3d_rms_m")-1)),
          ("CSMEngineIdleGuardPitch",engine_value("rpm800","guard","pitch_rms_deg")),
          ("CSMEngineIdleDeployedPitch",engine_value("rpm800","deployed","pitch_rms_deg")),
          ("CSMEngineRawRatio",engine_value("rpm2400","off","disp_3d_rms_m")/baseline),
          ("CSMEngineDeployedRatio",engine_value("rpm2400","deployed","disp_3d_rms_m")/baseline)):
        contents.append("\\newcommand{\\"+name+"}{"+f"{value:.3f}"+"}")
    for segment,name in (("low_start","LowStart"),("high","High"),("low_return","LowReturn")):
        value=df[(df.study=="roundtrip")&(df.method=="CSM")][f"seg_{segment}_disp_z_rms_m"].mean()
        contents.append("\\newcommand{\\CSMTransition"+name+"Heave}{"+f"{value:.3f}"+"}")
    (out/"csm-validation-tables-generated.tex-part").write_text("\n".join(contents)+"\n")
    return effects


def publish(results,out):
    out.mkdir(parents=True,exist_ok=True)
    protocol=json.loads((results/"protocol.json").read_text())
    if len(protocol["seed_triplets"])!=10:
        raise ValueError("publication tables require the full ten-seed cohort")
    df=pd.read_csv(results/"scores.csv")
    # The transition builder uses the released 8.5 m source, scaled to a 4 m
    # high endpoint. Preserve the source height separately from the scoring
    # reference height (also encoded in the simulator's H4 input filename).
    if "source_record_hs_m" not in df:
        df["source_record_hs_m"]=df.hs_m
    df.loc[df.study=="roundtrip","hs_m"]=4.0
    df.to_csv(results/"scores.csv",index=False,float_format="%.9g")
    if df.key.duplicated().any() or not (df.score_samples==180000).all():
        raise ValueError("duplicate rows or non-full-rate scoring window")
    finite=df[df.method.isin(KALMAN)][["disp_z_rms_m","disp_3d_rms_m","roll_rms_deg","pitch_rms_deg","yaw_rms_deg"]]
    if not np.isfinite(finite.to_numpy()).all():
        raise ValueError("nonfinite estimator results must be reported before publication")
    df["disp_3d_pct_hs"]=100*df.disp_3d_rms_m/df.hs_m
    e=engine_summary(df)
    e.to_csv(results/"engine_summary.csv",index=False,float_format="%.9g")
    effects=tables(df,e,out)
    nominal_figures(results,out)
    comparative_figure(df,out)
    engine_figures(e,out)
    transition_figure(results,out)
    summary={"replays":len(df),"paired_stationary_effects":effects,
             "transition_heights_m":{"low":1.5,"high":4.0,"released_high_source":8.5},
             "publication_source_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in (Path(__file__),ROOT/"plots/vector_export.py",ROOT/"plots/figure_decimation.py")},
             "files_sha256":{str(p.relative_to(results)):hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in [results/"protocol.json",results/"scores.csv",results/"engine_summary.csv",*sorted((results/"traces").glob("*.csv.gz"))]},
             "figures":{p.name:p.stat().st_size for p in sorted(out.glob("csm_validation_*.svg"))}}
    (results/"publication_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps({"replays":len(df),"paired_effects":effects,"figure_sizes":summary["figures"]},indent=2))
