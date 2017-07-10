#!/usr/bin/env python2.7
"""
Make some figures for the .tsv output of computeVariantsDistances.py
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string, math, copy
from collections import defaultdict
from Bio.Phylo.TreeConstruction import _DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')
import pylab
import networkx as nx
from collections import defaultdict
from toillib import robust_makedirs
from callVariants import alignment_sample_tag, alignment_region_tag, alignment_graph_tag, run
from callVariants import graph_path, sample_vg_path, g1k_vg_path, graph_path
from evaluateVariantCalls import defaultdict_set
from computeVariantsDistances import vcf_dist_header, read_tsv, write_tsv

# Set up the plot parameters
# Include both versions of the 1kg SNPs graph name
# copied from plotVariantComparison.sh

PLOT_PARAMS = [
    "--categories",
    "snp1kg",
    "snp1000g",
    "haplo1kg",
    "sbg",
    "cactus",
    "camel",
    "curoverse",
#    "debruijn-k31",
    "debruijn-k63",
    "level1",
    "level2",
    "level3",
    "prg",
    "refonly",
    "simons",
    "trivial",
    "vglr",
    "haplo1kg30",
    "haplo1kg50",
    "shifted1kg",
    "gatk3",
    "platypus",
    "g1kvcf",
    "freebayes",
    "samtools",
    "snp1kg_af001",
    "snp1kg_af010",
    "snp1kg_af100",
#    "snp1kg_kp",
    "haplo1kg30_af001",
    "haplo1kg30_af010",
    "haplo1kg30_af100",
    "haplo1kg50_af001",
    "haplo1kg50_af010",
    "haplo1kg50_af100",
#    "platinum",
    "freebayes_g",
    "snp1kg_norm",
    "snp1kg_plat",
    "--category_labels ",
    "1KG",
    "1KG",
    "\"1KG Haplo\"",
    "7BG",
    "Cactus",
    "Camel",
    "Curoverse",
 #   "\"De Bruijn 31\"",
    "\"De Bruijn 63\"",
    "Level1",
    "Level2",
    "Level3",
    "PRG",
    "Primary",
    "SGDP",
    "Unmerged",
    "VGLR",
    "\"1KG Haplo 30\"",
    "\"1KG Haplo 50\"",
    "Scrambled",
    "GATK3",
    "Platypus",
    "\"1000 Genomes\"",
    "FreeBayes",
    "Samtools",
    "\"1KG .001\"",
    "\"1KG .010\"",
    "\"1KG .100\"",
#    "\"1KG UF\"",
    "\"1KG Hap30 .001\"",
    "\"1KG Hap30 .010\"",
    "\"1KG Hap30 .100\"",
    "\"1KG Hap50 .001\"",
    "\"1KG Hap50 .010\"",
    "\"1KG Hap50 .100\"",
 #   "Platinum",
    "\"Freebayes VG\"",
    "\"1KG Norm\"",
    "\"1KG Plat\"",
    "--colors",
    "\"#fb9a99\"",
    "\"#fb9a99\"",
    "\"#fdbf6f\"",
    "\"#b15928\"",
    "\"#1f78b4\"",
    "\"#33a02c\"",
    "\"#a6cee3\"",
  #  "\"#e31a1c\"",
    "\"#ff7f00\"",
    "\"#FF0000\"",
    "\"#00FF00\"",
    "\"#0000FF\"",
    "\"#6a3d9a\"",
    "\"#000000\"",
    "\"#b2df8a\"",
    "\"#b1b300\"",
    "\"#cab2d6\"",
    "\"#00FF00\"",
    "\"#0000FF\"",
    "\"#FF0000\"",
    "\"#25BBD4\"",
    "\"#9E7C72\"",
    "\"#cab2d6\"",
    "\"#FF00FF\"",
    "\"#2F4F4F\"",
    "\"#C71585\"",
    "\"#663399\"",
    "\"#F4A460\"",
#    "\"#FA8072\"",
    "\"#556B2F\"",
    "\"#DEB887\"",
    "\"#800000\"",
    "\"#6A5ACD\"",
    "\"#C71585\"",
    "\"#FF6347\"",
#   "\"#119911\"",
    "\"#b1b300\"",
    "\"#fb4a44\"",
    "\"#fabf1f\"",
    "--font_size 20 --dpi 90"]

def make_filter_params(plot_params):
    # hack plot params to contain duplicate category
    p = []
    cat = None
    for param in plot_params:
        if param.startswith("--"):
            cat = param
            p.append(param)
        else:
            p.append(param)
            if cat == "--categories":
                p.append(param + ".filter")
            elif cat == "--category_labels ":
                p.append(".".format(param))
            else:
                p.append(param)
    num_cats = plot_params.index("--category_labels ") - plot_params.index("--categories") + 1
    p.append("--markers")
    for i in range(num_cats):
        p.append("o")
        p.append("+")
    return p

def name_map(options):
    """ make a dictionary from the list above """
    i = options.plot_params.index("--categories")
    j = options.plot_params.index("--category_labels ")

    names = dict()
    for k in range(i + 1, j):
        if options.plot_params[k][0:2] == "--":
            break
        names[options.plot_params[i + k]] = options.plot_params[j + k].replace("\"", "")

    # hack in sample-<name> and base-<name> maps
    keys = [k for k in names.keys()]
    for key in keys:
        names["base-{}".format(key)] = "Base {}".format(names[key])
        names["sample-{}".format(key)] = "Sample {}".format(names[key])
        
    return names

    
def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # General options
    parser.add_argument("comp_dir", type=str,
                        help="directory of comparison output written by computeVariantsDistances.py")
    parser.add_argument("--skip", type=str, default=None,
                        help="comma separated list of skip words to pass to plotHeatMap.py")

    parser.add_argument("--top", action="store_true",
                        help="print some zoom-ins too")
    parser.add_argument("--range", help="distance range to plot on either side of max f1 pr dot",
                        type=float, default=0.1)
    parser.add_argument("--max_steps", help="max points to display.  interpolate if more",
                        type=int, default=50)
    parser.add_argument("--totals", help="only do overall plots",
                        action="store_true")
    parser.add_argument("--name", default="Platinum Genomes", help="Name for comparison")                        
    parser.add_argument("--filter_comp_dir", type=str, default="None",
                        help="another comparison directory to show with comp_dir")

                            
    args = args[1:]

    return parser.parse_args(args)

def plot_kmer_comp(tsv_path, options):
    """ take a kmer compare table and make a 
    jaccard boxplot for the first column and a 
    recall / precision ploot for the 2nd and third column
    """
    out_dir = os.path.join(options.comp_dir, "comp_plots")
    robust_makedirs(out_dir)
    out_name = os.path.basename(os.path.splitext(tsv_path)[0])
    out_base_path = os.path.join(out_dir, out_name)
    sample = out_name.split("-")[-1].upper()
    region = out_name.split("-")[-2].upper()

    params = " ".join(options.plot_params)
    # jaccard boxplot
    jac_tsv = out_base_path + "_jac.tsv"
    awkstr = '''awk '{if (NR!=1) print $1 "\t" $2}' '''
    run("{} {} > {}".format(awkstr, tsv_path, jac_tsv))
    jac_png = out_base_path + "_jac.png"
    run("scripts/boxplot.py {} --save {} --title \"{} KMER Set Jaccard\" --x_label \"Graph\" --y_label \"Jaccard Index\" --x_sideways {}".format(jac_tsv, jac_png, region, params))

    # precision recall scatter plot
    acc_tsv = out_base_path + "_acc.tsv"
    awkstr = '''awk '{if (NR!=1) print $1 "\t" $4 "\t" $3}' '''
    run("{} {} > {}".format(awkstr, tsv_path, acc_tsv))
    acc_png = out_base_path + "_acc.png"
    run("scripts/scatter.py {} --save {} --title \"{} KMER Set Accuracy\" --x_label \"Recall\" --y_label \"Precision\" --width 12 --height 9 --lines {}".format(acc_tsv, acc_png, region, params))

def make_max_f1_tsv(acc_tsv_path, f1_tsv_path, f1_pr_tsv_path, f1_qual_tsv_path, f1_scat_tsv_path, options):
    """ flatten precision-recall tsv into single best f1 entry per graph """
    def f1(p, r):
        return 0 if p + r == 0 else 2. * ((p * r) / (p + r))
    max_f1 = defaultdict(int)
    max_pr = dict()
    max_qual = dict()
    with open(acc_tsv_path) as f:
        pr_file = [line for line in f]
    for i, line in enumerate(pr_file):
        toks = line.split()
        try:
            name, recall, precision, qual = toks[0], float(toks[1]), float(toks[2]), float(toks[3])
            max_f1[name] = max(f1(precision, recall), max_f1[name])
            if max_f1[name] == f1(precision, recall):
                max_pr[name] = (precision, recall, i)
                max_qual[name] = qual
        except:
            pass
    with open(f1_tsv_path, "w") as f1_file:
        for name, f1_score in max_f1.items():
            f1_file.write("{}\t{}\n".format(name, f1_score))
    with open(f1_pr_tsv_path, "w") as f1_pr_file, open(f1_scat_tsv_path, "w") as f1_scat_file:
        for name, pr_score in max_pr.items():
            best_f1_line = pr_file[pr_score[2]]
            best_recall, best_precision = float(best_f1_line.split()[1]), float(best_f1_line.split()[2])
            # write all points that are within range of max f1 on either both axes
            best_lines = []
            for i in range(0, len(pr_file)):
                line = pr_file[i]
                toks = line.split()
                if toks[0] == name:
                    recall, precision = float(toks[1]), float(toks[2])
                    if abs(recall - best_recall) <= options.range and abs(precision - best_precision) <= options.range:
                        best_lines.append((name, recall, precision))
                                    
            ds_step = max(1, int((len(best_lines)) / options.max_steps))
            for i, line in enumerate(best_lines):
                if i % ds_step == 0 or i == 0 or i == len(best_lines) - 1:                    
                    f1_pr_file.write("{}\t{}\t{}\n".format(line[0], line[1], line[2]))
                            
            # write the single max f1 point as precision recall
            toks = best_f1_line.split()
            recall, precision = float(toks[1]), float(toks[2])
            f1_scat_file.write("{}\t{}\t{}\n".format(name, recall, precision))

    with open(f1_qual_tsv_path, "w") as f1_qual_file:
        for name, qual_score in max_qual.items():
            f1_qual_file.write("{}\t{}\n".format(name, qual_score))
            

def plot_vcf_comp(tsv_path, options):
    """ take the big vcf compare table and make precision_recall plots for all the categories"""
    out_dir = os.path.join(options.comp_dir, "comp_plots")
    robust_makedirs(out_dir)
    out_name = os.path.basename(os.path.splitext(tsv_path)[0])
    sample = out_name.split("-")[-1].upper()
    region = out_name.split("-")[-2].upper()

    if options.totals is True and (sample.upper() != "COMBINED" or region.upper() != "TOTAL"):
        return
    
    def out_base_path(tag, label, extension):
        bd = tag if extension != ".tsv" else "tsv"
        ret = os.path.join(out_dir, bd, "-".join(out_name.split("-")[:-1]) + "-{}-{}-".format(sample, tag) + region) + "_" + label + extension
        robust_makedirs(os.path.dirname(ret))
        return ret

    params = " ".join(options.plot_params)

    # precision recall scatter plot
    header = vcf_dist_header(options)
    # strip qual
    header = header[:-1]
    for i in range(len(header) / 2):
        prec_idx = 2 * i
        rec_idx = prec_idx + 1
        qual_idx = len(header)
        print prec_idx, header[prec_idx], rec_idx, header[rec_idx]
        ptoks = header[prec_idx].split("-")
        rtoks = header[rec_idx].split("-")
        assert ptoks[1] == "Precision"
        assert rtoks[1] == "Recall"
        assert ptoks[:1] == rtoks[:1]
        comp_cat  = ptoks[0]
        if comp_cat not in ["TOT", "SNP", "INDEL"]:
            continue
        label = header[prec_idx].replace("Precision", "acc")
        acc_tsv = out_base_path("pr", label, ".tsv")
        print "Make {} tsv with cols {} {}".format(label, rec_idx, prec_idx)
        # +1 to convert to awk 1-base coordinates. +1 again since header doesnt include row_label col
        awkcmd = '''if (NR!=1) print $1 "\t" ${} "\t" ${} "\t" ${}'''.format(rec_idx + 2, prec_idx + 2, qual_idx + 2)
        awkstr = "awk \'{" + awkcmd + "}\'"
        run("{} {} > {}".format(awkstr, tsv_path, acc_tsv))
        acc_png = out_base_path("pr", label, ".png")
        title = ""
        if comp_cat == "TOT":
            title += "{} Accuracy".format(options.name)
        else:
            title += "{} {} Accuracy".format(options.name, comp_cat.title())
        if sample.upper() != "COMBINED" and "JOIN" not in sample.upper():
            title += ", " + sample.upper()
        if region.upper() != "TOTAL" and region.upper() != sample.upper():
            title += ", {}".format(region)
            
        cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 18 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --min_x -0.01 --max_x 1.01 --min_y -0.01 --max_y 1.01".format(acc_tsv, acc_png, title, params)
        print cmd
        os.system(cmd)

        #flatten to max f1 tsv and plot as bars
        f1_tsv = out_base_path("f1bar", label, ".tsv")
        f1_png = out_base_path("f1bar", label, ".png")
        f1_pr_tsv = out_base_path("f1pr", label, ".tsv")
        f1_pr_png = out_base_path("f1pr", label, ".png")
        f1_qual_tsv = out_base_path("f1qual", label, ".tsv")
        f1_qual_png = out_base_path("f1qual", label, ".png")
        f1_scat_tsv = out_base_path("f1scat", label, ".tsv")
        f1_scat_png = out_base_path("f1scat", label, ".png")


        make_max_f1_tsv(acc_tsv, f1_tsv, f1_pr_tsv, f1_qual_tsv, f1_scat_tsv, options)

        cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --annotate --no_legend".format(f1_pr_tsv, f1_pr_png, title, params)
        print cmd
        os.system(cmd)

        if not options.filter_comp_dir:
            cmd = "scripts/barchart.py {} --ascending --no_n --save {} --title \"{}\" --x_sideways --x_label \"Graph\" --y_label \"Max. F1-Score\" {}".format(f1_tsv, f1_png, title, params)
            print cmd
            os.system(cmd)
            cmd = "scripts/barchart.py {} --ascending --no_n --save {} --title \"{}\" --x_sideways --x_label \"Graph\" --y_label \"Quality for Max. F1-Score\" {}".format(f1_qual_tsv, f1_qual_png, title, params)
            print cmd
            os.system(cmd)
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --annotate --no_legend".format(f1_scat_tsv, f1_scat_png, title, params)
            print cmd
            os.system(cmd)

        if options.top is True:
            # top .25 f1pr scatter
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --min_x 0.746 --max_x 1.004 --min_y 0.746 --max_y 1.004 --annotate --no_legend".format(f1_pr_tsv, f1_pr_png.replace(".png", "_top25.png"), title, params)
            print cmd
            os.system(cmd)

            # top .50 f1pr scatter
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --min_x 0.496 --max_x 1.004 --min_y 0.496 --max_y 1.004 --annotate --no_legend".format(f1_pr_tsv, f1_pr_png.replace(".png", "_top50.png"), title, params)
            print cmd
            os.system(cmd)

            # top .65 f1pr scatter
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --min_x 0.646 --max_x 1.004 --min_y 0.646 --max_y 1.004 --annotate --no_legend".format(f1_pr_tsv, f1_pr_png.replace(".png", "_top65.png"), title, params)
            print cmd
            os.system(cmd)

            # top .70 f1pr scatter
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --min_x 0.696 --max_x 1.004 --min_y 0.696 --max_y 1.004 --annotate --no_legend".format(f1_pr_tsv, f1_pr_png.replace(".png", "_top70.png"), title, params)
            print cmd
            os.system(cmd)

            # top .85 f1pr scatter
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --min_x 0.846 --max_x 1.004 --min_y 0.846 --max_y 1.004 --annotate --no_legend".format(f1_pr_tsv, f1_pr_png.replace(".png", "_top85.png"), title, params)
            print cmd
            os.system(cmd)

            # top .90 f1pr scatter
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --min_x 0.896 --max_x 1.004 --min_y 0.896 --max_y 1.004 --annotate --no_legend".format(f1_pr_tsv, f1_pr_png.replace(".png", "_top90.png"), title, params)
            print cmd
            os.system(cmd)

            # top .95 f1pr scatter
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --min_x 0.946 --max_x 1.004 --min_y 0.946 --max_y 1.004 --annotate --no_legend".format(f1_pr_tsv, f1_pr_png.replace(".png", "_top95.png"), title, params)
            print cmd
            os.system(cmd)

        if not options.filter_comp_dir and options.top is True:
            # top 20
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 18 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --min_x 0.798 --max_x 1.002 --min_y 0.798 --max_y 1.002".format(acc_tsv, acc_png.replace(".png", "_top20.png"), title, params)
            print cmd
            os.system(cmd)
            # top 20
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 11 --height 5.5 {} --lines --no_n --line_width 1.5 --marker_size 5 --min_x 0.796 --max_x 1.004 --min_y 0.796 --max_y 1.004".format(acc_tsv, acc_png.replace(".png", "_top20_inset.png"), title, params)
            print cmd
            os.system(cmd)        
            # top 40
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 18 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --min_x 0.596 --max_x 1.004 --min_y 0.596 --max_y 1.004".format(acc_tsv, acc_png.replace(".png", "_top40.png"), title, params)
            print cmd
            os.system(cmd)
            # top .5 bar
            cmd = "scripts/barchart.py {} --ascending --no_n --save {} --title \"{}\" --x_sideways --x_label \"Graph\" --y_label \"Max. F1-Score\" {} --min 0.5".format(f1_tsv, f1_png.replace(".png", "_top50.png"), title, params)
            print cmd
            os.system(cmd)
            # top .6 bar
            cmd = "scripts/barchart.py {} --ascending --no_n --save {} --title \"{}\" --x_sideways --x_label \"Graph\" --y_label \"Max. F1-Score\" {} --min 0.6".format(f1_tsv, f1_png.replace(".png", "_top60.png"), title, params)
            print cmd
            os.system(cmd)
            # top .7 bar
            cmd = "scripts/barchart.py {} --ascending --no_n --save {} --title \"{}\" --x_sideways --x_label \"Graph\" --y_label \"Max. F1-Score\" {} --min 0.7".format(f1_tsv, f1_png.replace(".png", "_top70.png"), title, params)
            print cmd
            os.system(cmd)            
            # top .85 bar
            cmd = "scripts/barchart.py {} --ascending --no_n --save {} --title \"{}\" --x_sideways --x_label \"Graph\" --y_label \"Max. F1-Score\" {} --min 0.85".format(f1_tsv, f1_png.replace(".png", "_top85.png"), title, params)
            print cmd
            os.system(cmd)

            # top .95 bar
            cmd = "scripts/barchart.py {} --ascending --no_n --save {} --title \"{}\" --x_sideways --x_label \"Graph\" --y_label \"Max. F1-Score\" {} --min 0.95".format(f1_tsv, f1_png.replace(".png", "_top95.png"), title, params)
            print cmd
            os.system(cmd)
            

            # top .50 f1scat 
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --annotate --no_legend --min_x 0.496 --max_x 1.004 --min_y 0.496 --max_y 1.004".format(f1_scat_tsv, f1_scat_png.replace(".png", "_top50.png"), title, params)
            print cmd
            os.system(cmd)

            # top .60 f1scat 
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --annotate --no_legend --min_x 0.596 --max_x 1.004 --min_y 0.596 --max_y 1.004".format(f1_scat_tsv, f1_scat_png.replace(".png", "_top60.png"), title, params)
            print cmd
            os.system(cmd)

            # top .70 f1scat 
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --annotate --no_legend --min_x 0.696 --max_x 1.004 --min_y 0.696 --max_y 1.004".format(f1_scat_tsv, f1_scat_png.replace(".png", "_top70.png"), title, params)
            print cmd
            os.system(cmd)

            # top .75 f1scat 
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --annotate --no_legend --min_x 0.746 --max_x 1.004 --min_y 0.746 --max_y 1.004".format(f1_scat_tsv, f1_scat_png.replace(".png", "_top75.png"), title, params)
            print cmd
            os.system(cmd)

            # top .80 f1scat 
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --annotate --no_legend --min_x 0.796 --max_x 1.004 --min_y 0.796 --max_y 1.004".format(f1_scat_tsv, f1_scat_png.replace(".png", "_top85.png"), title, params)
            print cmd
            os.system(cmd)

        if options.top is True:
            
            # top .85 f1scat 
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --annotate --no_legend --min_x 0.846 --max_x 1.004 --min_y 0.846 --max_y 1.004".format(f1_scat_tsv, f1_scat_png.replace(".png", "_top85.png"), title, params)
            print cmd
            os.system(cmd)

            # top .90 f1scat 
            cmd = "scripts/scatter.py {} --save {} --title \"{}\" --x_label \"Recall\" --y_label \"Precision\" --width 9 --height 9 {} --lines --no_n --line_width 1.5 --marker_size 5 --annotate --no_legend --min_x 0.896 --max_x 1.004 --min_y 0.896 --max_y 1.004".format(f1_scat_tsv, f1_scat_png.replace(".png", "_top90.png"), title, params)
            print cmd
            os.system(cmd)            





def plot_heatmap(tsv, options):
    """ make a heatmap """
    out_dir = os.path.join(options.comp_dir, "heatmaps")
    robust_makedirs(out_dir)
    mat, col_names, row_names, row_label = read_tsv(tsv)
    names = name_map(options)

    for i in range(len(col_names)):
        if col_names[i] in names:
            col_names[i] = names[col_names[i]]
    for i in range(len(row_names)):
        if row_names[i] in names:
            row_names[i] = names[row_names[i]]

    if "_rename" in tsv:
        return
    fix_tsv = tsv.replace(".tsv", "_rename.tsv")
    write_tsv(fix_tsv, mat, col_names, row_names, row_label)

    out_hm = os.path.join(out_dir, os.path.basename(tsv).replace(".tsv", ".png"))
    ph_opts = "--skip {}".format(options.skip) if options.skip is not None else ""
    cmd = "scripts/plotHeatmap.py {} {} {}".format(fix_tsv, out_hm, ph_opts)
    print cmd
    os.system(cmd)

    cmd = "scripts/plotHeatmap.py {} {} {} --log_scale".format(fix_tsv, out_hm.replace(".png", "_log.png"), ph_opts)
    print cmd
    os.system(cmd)
    
    
def main(args):
    
    options = parse_args(args)

    # if we have a separated "filter" input directory, we double up all our options to include
    # .filter versions of everything.  
    options.plot_params = PLOT_PARAMS if not options.filter_comp_dir else make_filter_params(PLOT_PARAMS)

    # look through tsvs in comp_tables
    for tsv in glob.glob(os.path.join(options.comp_dir, "comp_tables", "*.tsv")):
        if "hm" in os.path.basename(tsv).split("-"):
            plot_heatmap(tsv, options)
        elif "kmer" in os.path.basename(tsv).split("-"):
            plot_kmer_comp(tsv, options)
        elif "vcf" in os.path.basename(tsv).split("-") or "sompy" in tsv.split("-") \
             or "happy" in tsv.split("-") or "vcfeval" in tsv.split("-"):
            if options.filter_comp_dir:
                filter_tsv = tsv.replace(options.comp_dir, options.filter_comp_dir)
                if os.path.isfile(filter_tsv):
                    # we merge on our "filter tsv" to the input tsv adding .filter to all the names
                    # in the first column
                    join_tsv = tsv.replace(".tsv", ".join.tsv")
                    shutil.copy2(tsv, join_tsv)
                    with open(join_tsv, "a") as jt, open(filter_tsv) as ft:
                        for line in ft:
                            toks = line.split()
                            if toks[0] != "Graph":
                                jt.write(toks[0]  + ".filter" + "\t" + "\t".join(toks[1:]) + "\n")
                    tsv = join_tsv
                
            plot_vcf_comp(tsv, options)
                                                                

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
