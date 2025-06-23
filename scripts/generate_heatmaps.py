#!/usr/bin/env python3
import os
import json
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, LinearSegmentedColormap
import matplotlib.patches as patches

def read_fasta(path):
    seq = ""
    with open(path) as f:
        for line in f:
            if line.startswith('>'): continue
            seq += line.strip()
    return seq

def read_msa(path):
    header = None
    seq = ""
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if header is None:
                    header = line[1:].strip()
                else:
                    break
            else:
                seq += line.strip()
    return header, seq

def find_alignment_region(msa_seq, target):
    N, tlen = len(msa_seq), len(target)
    for i in range(N):
        cnt, j = 0, i
        while j < N and cnt < tlen:
            if msa_seq[j] != '-': cnt += 1
            j += 1
        if cnt < tlen:
            break
        if msa_seq[i:j].replace('-', '') == target:
            return i, j - 1
    raise ValueError("Could not find concatenated seqA+seqB in MSA")

def build_maps(msa_seq, start, end, lenA):
    mapA, mapB = {}, {}
    cnt = 0
    for col in range(start, end+1):
        if msa_seq[col] != '-':
            cnt += 1
            if cnt <= lenA:
                mapA[col+1] = cnt
            else:
                mapB[col+1] = cnt - lenA
    return mapA, mapB

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--project', required=True,
                   help='Project name (folder under results/)')
    args = p.parse_args()
    proj = args.project

    # file paths
    paired_json = f"results/find_homologues/{proj}/paired_sequences.json"
    msa_file    = f"results/msas/{proj}/proteinAB.aln"
    coup_file   = f"results/coevolution/{proj}/coevolution_results.couplings"
    out_dir     = f"results/heatmaps/{proj}"
    os.makedirs(out_dir, exist_ok=True)

    # 1) load query sequences
    with open(paired_json) as f:
        entry = json.load(f)[0]
    if 'proteinA' in entry and 'proteinB' in entry:
        seqA = entry['proteinA']['seq']
        seqB = entry['proteinB']['seq']
        headerA = entry['proteinA'].get('id','queryAused')
        headerB = entry['proteinB'].get('id','queryBused')
    else:
        if isinstance(entry, (list,tuple)) and len(entry)>=2:
            seqA = read_fasta(entry[0])
            seqB = read_fasta(entry[1])
        elif isinstance(entry, dict) and 'A' in entry and 'B' in entry:
            seqA = read_fasta(entry['A'])
            seqB = read_fasta(entry['B'])
        else:
            raise ValueError(f"Cannot extract sequences from JSON entry: {entry}")
        headerA, headerB = 'queryAused','queryBused'
    lenA, lenB = len(seqA), len(seqB)
    # write used query FASTAs
    with open(f"{out_dir}/queryAused.fasta",'w') as fa:
        fa.write(f">{headerA}\n{seqA}\n")
    with open(f"{out_dir}/queryBused.fasta",'w') as fb:
        fb.write(f">{headerB}\n{seqB}\n")

    full_seq = seqA + seqB

    # 2) read MSA & find window
    _, msa_seq = read_msa(msa_file)
    start_col, end_col = find_alignment_region(msa_seq, full_seq)

    # 3) build maps
    mapA, mapB = build_maps(msa_seq, start_col, end_col, lenA)

    # 4) load couplings (i=col0,j=col2,value=col5)
    df = pd.read_csv(
        coup_file, sep=r'\s+', header=None,
        usecols=[0,2,5], names=['ref_i','ref_j','raw'],
        comment='#', dtype=str
    )
    df['ref_i']=pd.to_numeric(df['ref_i'],errors='coerce')
    df['ref_j']=pd.to_numeric(df['ref_j'],errors='coerce')
    df['raw']  =pd.to_numeric(df['raw'],errors='coerce')
    df.dropna(subset=['ref_i','ref_j','raw'], inplace=True)
    df['ref_i']=df['ref_i'].astype(int)
    df['ref_j']=df['ref_j'].astype(int)
    df['value']=df['raw'].abs()

    # precompute inter-protein list
    recs=[]
    for _,r in df.iterrows():
        i,j,v=r.ref_i,r.ref_j,r.value
        if i in mapA and j in mapB:
            recs.append((i,j,'A',mapA[i],'B',mapB[j],v))
        elif j in mapA and i in mapB:
            recs.append((i,j,'B',mapB[i],'A',mapA[j],v))
    mapped=pd.DataFrame(recs,columns=[
        'ref_i','ref_j','protein_i','pos_i','protein_j','pos_j','score'])
    mapped.to_csv(f"{out_dir}/couplings_mapped.csv",index=False)
    inter_df = mapped[['pos_i','pos_j','score']].rename(
        columns={'pos_i':'A_idx','pos_j':'B_idx'}
    )

    # --- 5) generate 10 heatmaps at top 5%,10%...50% of ALL coupling scores ---
    heat_pcts = list(range(5, 55, 5))  # [5,10,...50]
    total = lenA + lenB
    grey_blue = LinearSegmentedColormap.from_list('grey_blue',['lightgrey','blue'])
    for top_pct in heat_pcts:
        thr = np.percentile(df['value'], 100-top_pct)
        G = np.full((total,total), np.nan)
        for _,r in df.iterrows():
            i,j,v=r.ref_i,r.ref_j,r.value
            if   i in mapA: i0=mapA[i]-1
            elif i in mapB: i0=lenA+mapB[i]-1
            else: continue
            if   j in mapA: j0=mapA[j]-1
            elif j in mapB: j0=lenA+mapB[j]-1
            else: continue
            G[i0,j0]=v; G[j0,i0]=v
        G[G < thr] = np.nan
        m = np.ma.masked_invalid(G)
        pos = m.compressed()
        if pos.size == 0:
            continue
        vmin, vmax = pos.min(), pos.max()
        fig,ax = plt.subplots(figsize=(8,8))
        cax = ax.imshow(m,origin='lower',cmap=grey_blue,
                        norm=LogNorm(vmin=vmin,vmax=vmax))
        fig.colorbar(cax,ax=ax,label=f'Top {top_pct}% (log scale)')
        ax.set(xlabel='Residue idx (A+B)', ylabel='Residue idx (A+B)',
               title=f'Global heatmap (top {top_pct}%)')
        plt.tight_layout()
        fig.savefig(f"{out_dir}/global_heatmap_top{top_pct}.png",dpi=300)
        plt.close(fig)

    # --- 6) networks at fixed percentiles ---
    net_pcts = [95.0, 99.0, 99.9, 99.99]
    thr_vals = [np.percentile(df['value'], p) for p in net_pcts]
    def plot_network(thr,label,fname):
        ya, yb = np.linspace(0,1,lenA), np.linspace(0,1,lenB)
        pos_net = {f'A{i}':(0,ya[i-1]) for i in range(1,lenA+1)}
        pos_net.update({f'B{j}':(1,yb[j-1]) for j in range(1,lenB+1)})
        disN=(0,ya[130]) if lenA>=131 else None
        disC=(ya[676],ya[-1]) if lenA>=677 else None

        top = inter_df[inter_df['score']>=thr]
        A_regs = find_segs(top['A_idx'])
        B_regs = find_segs(top['B_idx'])

        fig,ax=plt.subplots(figsize=(10,6)); ax.axis('off')
        for x,y in pos_net.values(): ax.scatter(x,y,s=5,color='lightgrey',alpha=0.5)
        if disN:
            ax.add_patch(patches.Rectangle((-0.05,disN[0]),0.05,disN[1]-disN[0],
                                          color='lightgrey',alpha=0.3))
            ax.text(-0.1,np.mean(disN),'disordered N',va='center',ha='right',fontsize=8)
        if disC:
            ax.add_patch(patches.Rectangle((-0.05,disC[0]),0.05,disC[1]-disC[0],
                                          color='lightgrey',alpha=0.3))
            ax.text(-0.1,np.mean(disC),'disordered C',va='center',ha='right',fontsize=8)
        cmap=plt.cm.Reds
        norm=Normalize(vmin=thr,vmax=top['score'].max())
        for s,e in A_regs:
            ys=ya[s-1:e]
            ax.add_patch(patches.Rectangle((-0.01,ys.min()),0.01,ys.max()-ys.min(),
                                          color='grey',alpha=0.5))
            ax.text(-0.02,ys.mean(),f'{s}-{e}' if s!=e else str(s),
                    va='center',ha='right',fontsize=6)
        for s,e in B_regs:
            ys=yb[s-1:e]
            ax.add_patch(patches.Rectangle((1,ys.min()),0.01,ys.max()-ys.min(),
                                          color='grey',alpha=0.5))
            ax.text(1.02,ys.mean(),f'{s}-{e}' if s!=e else str(s),
                    va='center',ha='left',fontsize=6)
        for _,r in top.iterrows():
            ai,bi,v=int(r.A_idx),int(r.B_idx),r.score
            x1,y1=pos_net[f'A{ai}']; x2,y2=pos_net[f'B{bi}']
            ax.plot([x1,x2],[y1,y2],lw=1.5,color=cmap(norm(v)),alpha=0.8)
        for i in top['A_idx'].unique():
            x,y=pos_net[f'A{int(i)}']; ax.scatter(x,y,s=20,color='black')
        for j in top['B_idx'].unique():
            x,y=pos_net[f'B{int(j)}']; ax.scatter(x,y,s=20,color='black')
        ax.set_title(f'{label} ({len(top)} edges)',pad=10)
        plt.tight_layout()
        fig.savefig(f"{out_dir}/{fname}.png",dpi=300)
        plt.close(fig)

    # plot at fixed percentiles
    for pct,thr in zip(net_pcts,thr_vals):
        plot_network(thr,f'Top {100-pct:.3g}%',f'network_top{pct}')

    # --- 7) adaptive networks for ~20 and ~100 edges ---
    scores = np.sort(inter_df['score'].values)[::-1]
    def edge_thresh(n):
        if len(scores) >= n:
            return scores[n-1]
        return scores.min()
    thr20 = edge_thresh(20)
    thr100= edge_thresh(100)
    plot_network(thr20, f'~20 edges ({(inter_df["score"]>=thr20).sum()})', 'network_20')
    plot_network(thr100,f'~100 edges ({(inter_df["score"]>=thr100).sum()})','network_100')

    # --- 8) counts vs threshold sweep (1%->0.001%) ---
    def find_segs(idxs):  # reuse helper
        idxs=sorted(set(int(x) for x in idxs))
        if not idxs: return []
        segs=[]; s=idxs[0]; p=s
        for x in idxs[1:]:
            if x==p+1: p=x
            else: segs.append((s,p)); s=p=x
        segs.append((s,p)); return segs

    p_fracs = np.logspace(-2, -5, 15)  # 1%->0.001%
    percents = 100*(1-p_fracs)
    counts = [(inter_df['score'] >= np.percentile(df['value'], pct)).sum()
              for pct in percents]
    fig,ax=plt.subplots()
    ax.plot(p_fracs*100, counts,'o-')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('Top coupling percentile (%)')
    ax.set_ylabel('Inter-protein count')
    ax.set_title('Counts vs threshold')
    plt.tight_layout()
    fig.savefig(f"{out_dir}/counts_vs_threshold.png",dpi=300)
    plt.close(fig)

if __name__=='__main__':
    main()
