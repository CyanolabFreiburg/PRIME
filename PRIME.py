import argparse
import os.path
import sys
from Bio.SeqUtils import GC
import math
import random
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from scipy.stats import genextreme


# PRIME - Profile Motif Search based on TSS

def check_input_parameters(in_fasta, in_pssm, in_motif, in_tss, in_pssm_s1, in_pssm_s2, in_pseudo_count, in_num_bkg_runs, in_alpha, in_out_folder, in_distribution):
    message = ""
    if not os.path.isfile(in_fasta):
        message += "-> Please specify  --fasta  parameter!\n"
    if not os.path.isfile(in_tss):
        message += "-> Please specify  --tss  parameter!\n"
    if args.pssm == "":
        if not os.path.isfile(in_motif):
            message += "-> Please specify  --motif  parameter!\n"
    if args.motif == "":
        if not os.path.isfile(in_pssm):
            message += "-> Please specify  --pssm  parameter!\n"
    if os.path.isfile(in_pssm) and os.path.isfile(in_motif):
        message += "-> Please specify only one of the named parameter: --pssm  or  --motif!\n"
    if in_pseudo_count < 0:
        message += "-> Parameter  --pseudo_count  takes only values >= 0!\n"
    if in_num_bkg_runs <= 0:
        message += "-> Parameter  --num_bkg_runs  takes only values > 0!\n"
    if in_alpha <= 0 or in_alpha > 1:
        message += "-> Parameter  --alpha  takes only values between the range (0, 1] !\n"
    if not os.path.isdir(in_out_folder):
        message += "-> Parameter  --out_folder  is not defined!\n"
    if in_distribution != "norm" and in_distribution != "genextreme":
        message += "-> Parameter --distribution  should be   norm  or  genextreme"
    return message


def read_genome_file(in_genome_fasta):
    genome_str = "X"
    handle = open(in_genome_fasta, "r")
    for line in handle:
        line = line.rstrip()
        if not line.startswith(">"):
            genome_str += line.upper()
    handle.close()
    gc_content = GC(genome_str) / 100
    return genome_str, gc_content


def setup_pssm(pssm_file, motif_file, gc_content, in_pseudocounts):
    pssm_lookup = dict()

    if os.path.isfile(pssm_file):
        # read precomputed pssm from file
        handle = open(pssm_file, "r")
        for line in handle:
            line = line.rstrip()
            if not line.startswith("#"):
                line_arr = line.split()
                if not line_arr[0] in pssm_lookup:
                    pssm_lookup[line_arr[0]] = list()
                    pssm_lookup[line_arr[0]] = [float(x) for x in line_arr[1:]]
        handle.close()
    else:
        # build pssm from sequences
        # (1) build PFM
        pssm_lookup["A"] = list()
        pssm_lookup["C"] = list()
        pssm_lookup["G"] = list()
        pssm_lookup["T"] = list()
        init = "false"
        w_size = 0
        handle = open(motif_file, "r")
        for line in handle:
            line = (line.rstrip()).upper()
            # init arrays
            if init == "false":
                w_size = len(line)
                pssm_lookup["A"] = [0] * w_size
                pssm_lookup["C"] = [0] * w_size
                pssm_lookup["G"] = [0] * w_size
                pssm_lookup["T"] = [0] * w_size
                init = "true"
            pos = 0
            for letter in list(line.strip()):
                pssm_lookup[letter][pos] += 1
                pos += 1
        handle.close()

        # (2) build PPM ABS to RELATIVE values
        for pos in range(0, w_size):
            c_sum = 0
            for nuc in ["A", "C", "G", "T"]:
                # get column sum
                if pssm_lookup[nuc][pos] == 0:
                    pssm_lookup[nuc][pos] = in_pseudocounts
                else:
                    pssm_lookup[nuc][pos] += in_pseudocounts
                c_sum += pssm_lookup[nuc][pos]
            for nuc in ["A", "C", "G", "T"]:
                # get column sum + pseudocounts
                pssm_lookup[nuc][pos] = pssm_lookup[nuc][pos] / c_sum

        # (3) build final PSSM
        gc_single_con = gc_content / 2
        at_single_con = (1 - gc_content) / 2
        nucl_prob = {"A": at_single_con, "T": at_single_con, "C": gc_single_con, "G": gc_single_con}
        for pos in range(0, w_size):
            for nuc in ["A", "C", "G", "T"]:
                pssm_lookup[nuc][pos] = math.log2( (pssm_lookup[nuc][pos] / nucl_prob[nuc]) )
    return pssm_lookup


def compute_empirical_background_model_norm(in_genome, in_pssm, in_background_rand):
    # compute all possible pssm-scores
    random_genome_pssm_scores = list()
    pssm_len = len(in_pssm["A"])
    genome_size = len(in_genome)

    for i in range(0, in_background_rand):
        rand_start_pos = random.randrange(1, genome_size, 1)
        # compute PSSM Score from random selected start position
        temp_score = 0
        temp_pssm_pos = 0
        if rand_start_pos + pssm_len < genome_size:
            for i in range(rand_start_pos, rand_start_pos + pssm_len):
                if in_genome[i] in in_pssm:
                    temp_score += in_pssm[in_genome[i]][temp_pssm_pos]
                temp_pssm_pos += 1
            random_genome_pssm_scores.append(temp_score)
    # Fit a normal distribution to the data:
    mu, std = norm.fit(random_genome_pssm_scores)
    return mu, std


def compute_empirical_background_model_genextreme(in_genome, in_pssm, in_background_rand):
    # compute all possible pssm-scores
    random_genome_pssm_scores = list()
    pssm_len = len(in_pssm["A"])
    genome_size = len(in_genome)

    for i in range(0, in_background_rand):
        rand_start_pos = random.randrange(1, genome_size, 1)
        # compute PSSM Score from random selected start position
        temp_score = 0
        temp_pssm_pos = 0
        if rand_start_pos + pssm_len < genome_size:
            for i in range(rand_start_pos, rand_start_pos + pssm_len):
                if in_genome[i] in in_pssm:
                    temp_score += in_pssm[in_genome[i]][temp_pssm_pos]
                temp_pssm_pos += 1
            random_genome_pssm_scores.append(temp_score)
    # Fit a normal distribution to the data:
    c, loc, scale = genextreme.fit(random_genome_pssm_scores)
    return c, loc, scale, random_genome_pssm_scores


def fwd_helper_for_analyse_tss_file(in_tss_start_pos, in_pssm_s1, in_pssm_s2, in_pssm, in_genome):
    result = list()
    pssm_len = len(in_pssm["A"])

    if in_pssm_s1 > 0:
        motif_range = "[+" + str(in_pssm_s1) + ", "
    else:
        motif_range = "[" + str(in_pssm_s1) + ", "

    if in_pssm_s2 > 0:
        motif_range += "+" + str(in_pssm_s2) + "]"
    else:
        motif_range += str(in_pssm_s2) + "]"

    if in_pssm_s1 <= 0:
        rel_start_pos = in_tss_start_pos + in_pssm_s1 - pssm_len + 1
    else:
        rel_start_pos = in_tss_start_pos + in_pssm_s1
    if in_pssm_s2 <= 0:
        rel_end_pos = in_tss_start_pos + in_pssm_s2 - pssm_len + 2
    else:
        rel_end_pos = in_tss_start_pos + in_pssm_s2 + 1

    for i in range(rel_start_pos, rel_end_pos):
        tmp_start_pos = i
        tmp_end_pos = i + pssm_len - 1
        if tmp_start_pos >= 1:
            temp_score = 0
            temp_pssm_pos = 0
            temp_motif = ""
            for j in range(tmp_start_pos, tmp_end_pos + 1):
                if in_genome[j] in in_pssm:
                    temp_motif += str(in_genome[j])
                    temp_score += in_pssm[in_genome[j]][temp_pssm_pos]
                temp_pssm_pos += 1

            if in_pssm_s1 > 0:
                motif_dist = "+" + str(tmp_start_pos - in_tss_start_pos)
            else:
                motif_dist = "-" + str(in_tss_start_pos - tmp_end_pos)

            result.append([in_tss_start_pos, motif_range, motif_dist, temp_motif, tmp_start_pos, tmp_end_pos, "+", temp_score])
    return result


def rev_helper_for_analyse_tss_file(in_tss_start_pos, in_pssm_s1, in_pssm_s2, in_pssm, in_genome):
    result = list()
    pssm_len = len(in_pssm["A"])
    lookup = {"A": "T", "C": "G", "G": "C", "T": "A"}

    if in_pssm_s1 > 0:
        motif_range = "[+" + str(in_pssm_s1) + ", "
    else:
        motif_range = "[" + str(in_pssm_s1) + ", "

    if in_pssm_s2 > 0:
        motif_range += "+" + str(in_pssm_s2) + "]"
    else:
        motif_range += str(in_pssm_s2) + "]"

    if in_pssm_s1 <= 0:
        rel_end_pos = in_tss_start_pos + (-1 * in_pssm_s1) + 1
    else:
        rel_end_pos = in_tss_start_pos - in_pssm_s1 - pssm_len + 2
    if in_pssm_s2 <= 0:
        rel_start_pos = in_tss_start_pos + (-1 * in_pssm_s2)
    else:
        rel_start_pos = in_tss_start_pos - in_pssm_s2 - pssm_len + 1

    for i in range(rel_start_pos, rel_end_pos):
        tmp_start_pos = i
        tmp_end_pos = i + pssm_len - 1

        if tmp_end_pos <= len(in_genome):
            temp_score = 0
            temp_pssm_pos = pssm_len - 1
            temp_motif = ""
            for j in range(tmp_start_pos, tmp_end_pos + 1):
                if lookup[in_genome[j]] in in_pssm:
                    compl_nucl = lookup[in_genome[j]]
                    temp_motif += str(compl_nucl)
                    # Calculate PSSM score backwards!
                    temp_score += in_pssm[compl_nucl][temp_pssm_pos]
                temp_pssm_pos -= 1
            temp_motif = temp_motif[::-1]

            if in_pssm_s1 > 0:
                motif_dist = "+" + str(in_tss_start_pos - tmp_start_pos)
            else:
                motif_dist = "-" + str(tmp_end_pos - in_tss_start_pos)

            result.append([in_tss_start_pos, motif_range, motif_dist, temp_motif, tmp_start_pos, tmp_end_pos, "-", temp_score])
    return result


def analyse_tss_file(in_tss_args, in_pssm_s1_args, in_pssm_s2_args, in_pssm, in_genome):
    master_table = list()
    handle = open(in_tss_args, "r")
    count = 0
    ensure_unique_tss = dict()

    # autodetection to get the correct order of the input parameters -x and -y
    if int(in_pssm_s1_args) >= int(in_pssm_s2_args):
        swap = in_pssm_s2_args
        in_pssm_s2_args = in_pssm_s1_args
        in_pssm_s1_args = swap

    for line in handle:
        line = line.rstrip()
        if count > 0:
            line_arr = line.split()
            # Auto detect input format => customized table
            if len(line_arr) <= 10:
                tss_start_pos = int(line_arr[0])
                orientation = line_arr[1].upper()
            # Auto detect input format => ReadXplorer format
            if len(line_arr) > 10:
                tss_start_pos = int(line_arr[0])
                orientation = line_arr[3].upper()
            if not tss_start_pos in ensure_unique_tss:
                ensure_unique_tss[tss_start_pos] = 0
                if orientation == "FWD" or orientation == "+":
                    master_table.append(fwd_helper_for_analyse_tss_file(tss_start_pos, in_pssm_s1_args, in_pssm_s2_args, in_pssm, in_genome))
                else:
                    master_table.append(rev_helper_for_analyse_tss_file(tss_start_pos, in_pssm_s1_args, in_pssm_s2_args, in_pssm, in_genome))
        count += 1
    handle.close()
    final_master_table = list()
    # reformat master table 3d -> 2d table
    for entry in master_table:
        for sub_entry in entry:
            final_master_table.append(sub_entry)
    return final_master_table


def calculate_adjusted_p_values_norm(in_master_table, in_mu, in_std, in_alpha):
    raw_p_values = list()
    for i in range(0, len(in_master_table)):
        tmp_p_val = norm(in_mu, in_std).sf(in_master_table[i][7])
        raw_p_values.append(tmp_p_val)
        master_table[i].append(tmp_p_val)
    # adjust p-values
    if len(raw_p_values) >= 2:
        adjusted_p_values = multipletests(raw_p_values, alpha=in_alpha, method='fdr_bh', is_sorted=False)
        for i in range(0, len(adjusted_p_values[1])):
            in_master_table[i].append(adjusted_p_values[1][i])
    else:
        for i in range(0, len(in_master_table)):
            in_master_table[i].append("na")
    return in_master_table


def calculate_adjusted_p_values_genextreme(in_master_table, c, loc, scale, in_alpha):
    raw_p_values = list()
    for i in range(0, len(in_master_table)):
        tmp_p_val = genextreme.sf(in_master_table[i][7], c, loc, scale)
        raw_p_values.append(tmp_p_val)
        master_table[i].append(tmp_p_val)
    # adjust p-values
    if len(raw_p_values) >= 2:
        adjusted_p_values = multipletests(raw_p_values, alpha=in_alpha, method='fdr_bh', is_sorted=False)
        for i in range(0, len(adjusted_p_values[1])):
            in_master_table[i].append(adjusted_p_values[1][i])
    else:
        for i in range(0, len(in_master_table)):
            in_master_table[i].append("na")
    return in_master_table


def build_final_tables(in_out_folder, in_master_table, in_alpha, in_description, in_pssm, in_gc, in_ps_counts, in_tss_list):
    # write csv file
    if in_out_folder.endswith('/'):
        path_csv = in_out_folder + str("motif-master-table.csv")
        path_gff = in_out_folder + str("motif.gff")
        path_pssm = in_out_folder + str("pssm.txt")
        path_tss = in_out_folder + str("tss.gff")
    else:
        path_csv = in_out_folder + "/" + str("motif-master-table.csv")
        path_gff = in_out_folder + "/" + str("motif.gff")
        path_pssm = in_out_folder + "/" + str("pssm.txt")
        path_tss = in_out_folder + "/" + str("tss.gff")

    handle = open(path_csv, "w")
    # write header line
    header = "#TSS-POS" + "\t" + "MOTIF-RANGE" + "\t" + "MOTIF-DISTANCE-TO-TSS" + "\t" + "MOTIF" + "\t" + "START-POS" + "\t" + "END-POS" + \
             "\t" + "STRAND" + "\t" + "PSSM-SCORE" + "\t" + "RAW_P-VALUE" + "\t" + "ADJUSTED_P-VALUE" + "\n"
    handle.write(header)
    # write data
    for entry in in_master_table:
        out_str = str("\t".join(map(str, entry))) + "\n"
        handle.write(out_str)
    handle.close()

    # write gff file and perform adjusted p-value dependent coloring
    handle = open(path_gff, "w")
    handle.write("##gff-version 3\n")
    for entry in in_master_table:
        if entry[9] < in_alpha:
            color = "255 10 225"
            last_col = "colour=" + str(color) + ";" + "note=" + "PSSM-Score " + str(round(entry[7], 6)) + " raw_p-Value " \
                       + str(round(entry[8], 6)) + " adjusted_p-Value " + str(round(entry[9], 6))
            out_str = "PRIME" + "\t" + "PRIME" + "\t" + str(in_description) + "\t" + str(entry[4]) + "\t" + \
                      str(entry[5]) + "\t" + "." + "\t" + str(entry[6]) + "\t" + "." + "\t" + str(last_col) + "\n"
            handle.write(out_str)
    handle.close()

    # write PSSM model to file
    handle = open(path_pssm, "w")
    # header line
    header = "# GC-Content: " + str(round(in_gc * 100, 2)) + "%  Used Pseudocount Value: " + str(in_ps_counts) + "\n"
    handle.write(header)
    letter_order = ["A", "C", "G", "T"]
    for letter in letter_order:
        entry = str(letter) + "\t" + "\t".join(map(str, in_pssm[letter])) + "\n"
        handle.write(entry)
    handle.close()

    # write TSS-File in GFF-Format
    # todo use the internal build up master table instead of reading the input file again
    handle_out = open(path_tss, "w")
    handle_in = open(in_tss_list, "r")
    count = 0
    unique_positions = dict()
    for line in handle_in:
        line = line.rstrip()
        line_arr = line.split()
        if count > 0 and not int(line_arr[0]) in unique_positions:
            orientation = line_arr[3].upper()
            if orientation == "FWD":
                orientation = "+"
            else:
                orientation = "-"
            out_str = "TSS-Start" + "\t" + "TSS-Start" + "\t" + "TSS" + "\t" + str(line_arr[0]) + "\t" + str(line_arr[0]) + \
                      "\t" + "." + "\t" + str(orientation) + "\t" + "." + "\t" + "colour=7 36 216;" + "\n"
            handle_out.write(out_str)
            unique_positions[int(line_arr[0])] = 0
        else:
            count += 1
    handle_out.close()
    handle_in.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="GENOME file in FASTA format <FILE>", type=str, default="")
    parser.add_argument("-p", "--pssm", help="PSSM file, otherwise it will be computed <FILE> [optional]", type=str, default="")
    parser.add_argument("-m", "--motif", help="Hand selected motif patterns <FILE> [optional if --pssm is used, otherwise it is mandatory!", type=str, default="")
    parser.add_argument("-s", "--tss", help="ReadXplorer Transcriptional Start Sites <FILE> or customized input table, "
                                            "first col. TSS-Start position, second col. Orientation and col. "
                                            "three to ten are optional for any typ of information", type=str, default="")
    parser.add_argument("-x", "--pssm_s1", help="PSSM: Distance 1 to related TSS", type=int, default=-15)
    parser.add_argument("-y", "--pssm_s2", help="PSSM: Distance 2 to detected TSS", type=int, default=-5)
    parser.add_argument("-c", "--pseudo_count", help="Value for pseudocounts. Value is added to zero as well as non-zero values", type=float, default=0.7)
    parser.add_argument("-n", "--num_bkg_runs", help="Number of randomly picked PSSM scores from genome to create a suitable background model", type=int, default=100000)
    parser.add_argument("-a", "--alpha", help="Alpha level used for the motif.gff file [0.0-1.0]", type=float, default=0.05)
    parser.add_argument("-o", "--out_folder", help="Folder for the final results", type=str, default="")
    parser.add_argument("-d", "--description", help="Feature type; only used for the GFF file", type=str, default="motif")
    parser.add_argument("-k", "--distribution", help="PSSM values are normal distributed or follows a generalized extreme value distribution; [norm, genextreme]", type=str, default="norm")
    args = parser.parse_args()

    # *** PREPROCESSING ***
    # check input for plausibilty
    in_message = check_input_parameters(args.fasta, args.pssm, args.motif, args.tss, args.pssm_s1, args.pssm_s2,
                                        args.pseudo_count, args.num_bkg_runs, args.alpha, args.out_folder, args.distribution)
    if in_message != "":
        sys.stderr.write(str(in_message) + "\n")
        exit()

    # *** MAIN PROCESSING STEPS ***
    # read GENOME FASTA file -> genome goes from 1 to n!
    genome_str, gc_content = read_genome_file(args.fasta)

    # read or build PSSM
    pssm = setup_pssm(args.pssm, args.motif, gc_content, args.pseudo_count)

    # compute background model
    if args.distribution == "norm":
        mu, std = compute_empirical_background_model_norm(genome_str, pssm, args.num_bkg_runs)
    else:
        c, loc, scale, data_b = compute_empirical_background_model_genextreme(genome_str, pssm, args.num_bkg_runs)

    # read TSS File and compute PSSM scores
    master_table = analyse_tss_file(args.tss, args.pssm_s1, args.pssm_s2, pssm, genome_str)

    # calculate p-values from PSSM scores
    if args.distribution == "norm":
        master_table = calculate_adjusted_p_values_norm(master_table, mu, std, args.alpha)
    else:
        master_table = calculate_adjusted_p_values_genextreme(master_table, c, loc, scale, args.alpha)

    #rv = genextreme(c=c, loc=loc, scale=scale)
    #x = np.linspace(-20, 20, len(data_b))  # Generate 50 values from 0-40
    #y_pdf = rv.pdf(x)
    #plt.plot(x, y_pdf, 'ro-', label='PDF')
    #plt.hist(data_b, bins=50, density=True)
    #plt.legend()
    #plt.show()

    # *** PROCESSING DATA FOR OUTPUT ***
    # transfer master table in csv + gff-format
    build_final_tables(args.out_folder, master_table, args.alpha, args.description, pssm, gc_content, args.pseudo_count, args.tss)