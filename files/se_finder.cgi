#!/usr/bin/python

#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License.

import cgi

import numpy

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def bases_to_num(base):

    lbase = base.lower()
    b_up_num, b_low_num = list(), list()
    flag = 0
    
    for i in lbase:
        flag = 1
        if i == 'a':
            b_up_num.append(0)
            b_low_num.append(1)
        elif i == 't':
            b_up_num.append(1)
            b_low_num.append(0)
        elif i == 'c':
            b_up_num.append(2)
            b_low_num.append(3)
        elif i == 'g':
            b_up_num.append(3)
            b_low_num.append(2)
        else:
            flag = 0
            break
        
    b_low_num.reverse()
    return b_up_num, b_low_num, flag

def num_to_base(num):
    base = ''
    for i in num:
        if i == 0:
            base = base + 'a'
        elif i == 1:
            base = base + 't'
        elif i == 2:
            base = base + 'c'
        elif i == 3:
            base = base + 'g'

    return base
            


def h_bond(b_up_num, b_low_num):

    HSCORE = [[2.28, 0.69, 2.64, 1.89, 0.61, -1.00, 0.88, 0.14, 2.64, 1.33, 3.01, 2.26, 1.89, 0.74, 2.26, 1.50],
              [0.61, 2.23, 2.52, 1.63, -0.88, 0.69, 0.73, 0.07, 0.77, 2.59, 2.88, 2.00, 0.02, 1.84, 2.13, 1.25],
              [2.64, 2.52, 2.77, 0.17, 0.77, 0.64, 1.33, -1.44, 3.01, 2.88, 3.14, 0.47, 2.26, 2.13, 2.38, -0.52],
              [1.89, 1.74, 0.43, 1.58, 0.02, 0.71, -1.28, -0.13, 2.26, 2.10, 0.79, 1.94, 1.51, 1.35, 0.11, 1.19],
              [0.69, -0.58, 0.92, 0.42, 2.23, 0.68, 2.59, 1.84, 2.52, 0.97, 2.88, 2.13, 1.74, 0.43, 2.10, 1.35],
              [-1.00, 0.68, 0.75, 0.34, 0.69, 2.18, 2.47, 1.58, 0.64, 2.47, 2.76, 1.87, 0.71, 1.69, 1.98, 0.74],
              [1.33, 0.97, 1.05, -1.30, 2.59, 2.47, 2.72, 0.45, 2.88, 2.75, 3.01, 0.62, 2.10, 1.97, 2.23, 0.08],
              [0.74, 0.43, -1.45, 0.44, 1.84, 1.69, -0.12, 1.53, 2.13, 1.97, 0.62, 1.81, 1.35, 0.52, -0.47, 1.03],
              [2.64, 0.92, 3.01, 2.26, 2.52, 0.75, 2.88, 2.13, 2.77, 1.05, 3.14, 2.38, 0.43, -1.45, 0.75, 0.03],
              [0.88, 2.59, 2.88, 2.00, 0.73, 2.47, 2.76, 1.87, 1.33, 2.72, 3.01, 2.13, -1.28, -0.12, 0.40, -0.32],
              [3.01, 2.88, 3.14, 0.81, 2.88, 2.76, 3.01, 0.98, 3.14, 3.01, 3.26, 0.79, 0.79, 0.62, 0.70, -1.84],
              [2.26, 2.10, 0.75, 1.94, 2.13, 1.98, 0.40, 1.82, 2.38, 2.23, 0.70, 2.07, 0.11, -0.47, -2.17, -0.11],
              [1.89, 0.42, 2.26, 1.50, 1.63, 0.34, 2.00, 1.24, 0.17, -1.30, 0.81, -0.25, 1.58, 0.44, 1.94, 1.19],
              [0.14, 1.84, 2.13, 1.24, 0.07, 1.58, 1.87, 1.15, -1.44, 0.45, 0.98, -0.59, -0.13, 1.53, 1.82, 0.93],
              [2.26, 2.13, 2.38, -0.25, 2.00, 1.87, 2.13, -0.59, 0.47, 0.62, 0.79, -2.24, 1.94, 1.81, 2.07, -1.11],
              [1.50, 1.35, 0.03, 1.19, 1.25, 0.74, -0.32, 0.93, -0.52, 0.08, -1.84, -1.11, 1.19, 1.03, -0.11, 0.87]]

    delG_hbond = 0
    n = min([len(b_up_num), len(b_low_num)])   #have to warn user to calculating mismatch with uneven length
    for i in range(n-1):
        a = b_up_num[i]
        b = b_up_num[i+1]
        c = b_low_num[n-i-1]
        d = b_low_num[n-i-2]
        delG_hbond = delG_hbond + HSCORE[a*4+b][c*4+d]

    return delG_hbond
        
def stack(opt, left1, left2, right1, right2):
    if opt == "off":
        STACKSCORE = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    elif opt == "forian":
        STACKSCORE = [[-0.8, -0.5, -0.3, -0.9],
                      [-0.5, -0.4, -0.2, -0.5],
                      [-0.3, -0.2, -0.6, -0.6],
                      [-0.9, -0.5, -0.6, -1.1]]
    elif opt == "protozanova":
        STACKSCORE = [[-1.11, -1.34, -1.81, -1.06],   
                      [-0.19, -1.11, -1.43, -0.55],
                      [-0.55, -1.06, -1.44, -0.91],
                      [-1.43, -1.81, -2.17, -1.44]]

    return STACKSCORE[left1][left2] + STACKSCORE[right1][right2]
    

def delG_cal(se_up, se_low, se_left, se_right, se_type, stack_opt, delG_lig):
    se_left_up, se_left_low, flag = bases_to_num(se_left)
    se_right_up, se_right_low, flag = bases_to_num(se_right)
    if se_type == "5":
        delG_se = h_bond(se_up, se_low) + stack(stack_opt,
                        se_left_up[0], se_up[0], se_right_low[0], se_low[0]) + delG_lig
    elif se_type == "3":
        delG_se = h_bond(se_up, se_low) + stack(stack_opt,
                        se_low[-1], se_left_low[0], se_up[-1], se_right_up[0]) + delG_lig

    return delG_se


def delG_mm_cal(se_up, se_low, se_left, se_right, se_type, stack_opt, delG_lig):
    se_left_up, se_left_low, flag = bases_to_num(se_left)
    se_right_up, se_right_low, flag = bases_to_num(se_right)
    if se_type == "5":
        uu = h_bond(se_up, se_up) + stack(stack_opt, se_right_low[0], se_up[0],
                    se_right_low[0], se_up[0]) + delG_lig
        ll = h_bond(se_low, se_low) + stack(stack_opt, se_left_up[0], se_low[0],
                    se_left_up[0], se_low[0]) + delG_lig
    elif se_type == "3":
        uu = h_bond(se_up, se_up) + stack(stack_opt, se_up[-1], se_left_low[0],
                    se_up[-1], se_left_low[-1]) + delG_lig
        ll = h_bond(se_low, se_low) + stack(stack_opt, se_low[-1], se_right_up[0],
                    se_low[-1], se_right_up[0]) + delG_lig
    uu_gap = h_bond(se_up[0:-1], se_up[0:-1]) + delG_lig
    ll_gap = h_bond(se_low[0:-1], se_low[0:-1]) + delG_lig
    ul_gap = h_bond(se_up[0:-1], se_low[0:-1]) + delG_lig
    return tuple([uu, ll, uu_gap, ll_gap, ul_gap])

def score_cal(delG, delG_mm):
    mm_mean = numpy.mean(delG_mm)
    mm_std = numpy.std(delG_mm)
    score = ((-1)*delG - (-1)*mm_mean)/mm_std
    '''
    score = 0
    score_each = numpy.zeros(len(delG_mm))
    for i in range(len(delG_mm)):
        score_each[i] = ((-1)*delG - (-1)*delG_mm[i])/((-1)*delG)
        if (i == 0 or i == 1):
            score = score + score_each[i]
        else:
            score = score + 0.2*score_each[i]'''
    return score

def score_cal_2se(dna1_delG, dna1_delG_mm, dna2_delG, dna2_delG_mm, dna12_delG_mm):
    dna1_mm_mean = numpy.mean(numpy.concatenate((dna1_delG_mm, dna12_delG_mm)))
    dna1_mm_std = numpy.std(numpy.concatenate((dna1_delG_mm, dna12_delG_mm)))
    dna2_mm_mean = numpy.mean(numpy.concatenate((dna2_delG_mm, dna12_delG_mm)))
    dna2_mm_std = numpy.std(numpy.concatenate((dna2_delG_mm, dna12_delG_mm)))
    score = ((-1)*dna1_delG - (-1)*dna1_mm_mean)/dna1_mm_std + ((-1)*dna2_delG - (-1)*dna2_mm_mean)/dna2_mm_std
                            
    '''
    score = 0
    delG = (dna1_delG + dna2_delG)/2
    score_each = numpy.zeros(len(delG_mm))
    for i in range(len(delG_mm)):
        score_each[i] = ((-1)*delG - (-1)*delG_mm[i])/((-1)*delG)
        if (i < 4):
            score = score + score_each[i]
        else:
            score = score + 0.2*score_each[i]'''
    return score
    
def delG_2se_mm_cal(se1_up, se1_low, se2_up, se2_low, se1_left, se1_right,
                    se2_left, se2_right, se_type, stack_opt, delG_lig):
    se1_left_up, se1_left_low, flag = bases_to_num(se1_left)
    se1_right_up, se1_right_low, flag = bases_to_num(se1_right)
    se2_left_up, se2_left_low, flag = bases_to_num(se2_left)
    se2_right_up, se2_right_low, flag = bases_to_num(se2_right)
    if se_type == "5":
        u1u2 = h_bond(se1_up, se2_up) + stack(stack_opt, se2_right_low[0], se1_up[0],
                      se1_right_low[0], se2_up[0]) + delG_lig
        l1l2 = h_bond(se1_low, se2_low) + stack(stack_opt, se2_left_up[0], se1_low[0],
                      se1_left_up[0], se2_low[0]) + delG_lig
        u1l2 = h_bond(se1_up, se2_low) + stack(stack_opt, se2_left_up[0], se1_up[0],
                      se1_right_low[0], se2_low[0]) + delG_lig
        l1u2 = h_bond(se1_low, se2_up) + stack(stack_opt, se1_left_up[0], se2_up[0],
                      se2_right_low[0], se1_low[0]) + delG_lig
    elif se_type == "3":
        u1u2 = h_bond(se1_up, se2_up) + stack(stack_opt, se1_up[-1], se2_left_low[0],
                      se2_up[-1], se1_left_low[0]) + delG_lig
        l1l2 = h_bond(se1_low, se2_low) + stack(stack_opt, se1_low[-1], se2_right_up[0],
                      se2_low[-1], se1_right_up[0]) + delG_lig
        u1l2 = h_bond(se1_up, se2_low) + stack(stack_opt, se1_up[-1], se2_right_up[0],
                      se2_low[-1], se1_left_low[0]) + delG_lig
        l1u2 = h_bond(se1_low, se2_up) + stack(stack_opt, se2_up[-1], se1_right_up[0],
                      se1_low[-1], se2_left_low[0]) + delG_lig
    u1u2_gap = h_bond(se1_up[0:-1], se2_up[0:-1]) + delG_lig
    l1l2_gap = h_bond(se1_low[0:-1], se2_low[0:-1]) + delG_lig
    u1l2_gap = h_bond(se1_up[0:-1], se2_low[0:-1]) + delG_lig
    l1u2_gap = h_bond(se1_low[0:-1], se2_up[0:-1]) + delG_lig
    return tuple([u1u2, l1l2, u1l2, l1u2, u1u2_gap, l1l2_gap, u1l2_gap, l1u2_gap])

                      
def error_output(flag):
    print "Error:"
    if flag == 0:
        print "DNA sequence is missing."
    elif flag == 1:
        print "Please only enter A/C/T/G for DNA sequence 1."
    elif flag == 2:
        print "Please only enter A/C/T/G for sticky end 2."
    elif flag == 3:
        print "Please select 5' or 3' overhang type."
    elif flag == 4:
        print "Invalid custom free energy for the repair of phosphodiester bonds."
    elif flag == 5:
        print "DNA sequence is too short."


def good_output(opt, se_score, display_row):
    if opt == 0:
        dna1_se, dna1_score = se_score            
    elif opt == 1:
        dna1_se, dna1_score, dna2_se, dna2_score = se_score
    elif opt == 2:
        dna1_se, dna1_score, dna2_se, dna2_score, dna12_score = se_score

    dna1_table = numpy.zeros([len(dna1_se), 3+len(dna1_se[0,:])])
    for i in range(len(dna1_se)):
        dna1_table[i,0] = i + 2
        dna1_table[i,1] = i + len(dna1_se[i])+1
        dna1_table[i,2:2+len(dna1_se[0,:])] = dna1_se[i]
        dna1_table[i,-1] = dna1_score[i]
    dna1_sorted = dna1_table[(-dna1_table[:,-1]).argsort()]
    if opt > 0:
        dna2_table = numpy.zeros([len(dna2_se), 3+len(dna2_se[0,:])])
        for i in range(len(dna2_se)):
            dna2_table[i,0] = i + 2
            dna2_table[i,1] = i + len(dna2_se[i]) + 1
            dna2_table[i,2:2+len(dna2_se[0,:])] = dna2_se[i]
            dna2_table[i,-1] = dna2_score[i]
        dna2_sorted = dna2_table[(-dna2_table[:,-1]).argsort()]

        if opt == 2:
            dna12_table = numpy.zeros([len(dna12_score), 5 + len(dna1_se[0,:]) + len(dna2_se[0,:])])
            for i in range(len(dna1_score)):
                for j in range(len(dna2_score)):
                    dna12_table[i*len(dna2_score)+j,0] = i+2
                    dna12_table[i*len(dna2_score)+j,1] = i+len(dna1_se[i])+1
                    dna12_table[i*len(dna2_score)+j,2:2+len(dna1_se[0,:])] = dna1_se[i]
                    dna12_table[i*len(dna2_score)+j,2+len(dna1_se[0,:])] = j+2
                    dna12_table[i*len(dna2_score)+j,3+len(dna1_se[0,:])] = j+len(dna2_se[j])+1
                    dna12_table[i*len(dna2_score)+j,4+len(dna1_se[0,:]):4+len(dna1_se[0,:])+len(dna2_se[0,:])] = dna2_se[j]
                    dna12_table[i*len(dna2_score)+j,-1] = dna12_score[i*len(dna2_score) + j]
            dna12_sorted = dna12_table[(-dna12_table[:,-1]).argsort()]
    
    print """
        <h3>Sticky ends from DNA sequence 1:</h3>
        <table border="1">
        <tr>
        <th>Base number</th>
        <th>Sticky end</th>
        <th>Score</th>
        </tr>
    """
    for i in range(min(display_row,len(dna1_se))):
        print """
            <tr>
            <td>%.0f to %.0f</td>
            <td>%s</td>
            <td>%.2f</td>
            </tr>
        """ %(dna1_sorted[i,0], dna1_sorted[i,1], num_to_base(dna1_sorted[i,2:2+len(dna1_se[0,:])]), dna1_sorted[i,-1])
    print '</table><br/>'

    if opt > 0:
        print """
            <h3>Sticky ends from DNA sequence 2:</h3>
            <table border="1">
            <tr>
            <th>Base number</th>
            <th>Sticky end</th>
            <th>Score</th>
            </tr>
        """
        for i in range(min(display_row,len(dna2_se))):
            print """
                <tr>
                <td>%.0f to %.0f</td>
                <td>%s</td>
                <td>%.2f</td>
                </tr>
            """ %(dna2_sorted[i,0], dna2_sorted[i,1], num_to_base(dna2_sorted[i,2:2+len(dna2_se[0,:])]), dna2_sorted[i,-1])
        print '</table><br/>'
        if opt == 2:
            print """
                <h3>Combined score for sticky ends from DNA sequence 1 and 2:</h3>
                <table border="1">
                <tr>
                <th>Base number 2</th>
                <th>Sticky end 1</th>
                <th>Base number 2</th>
                <th>Sticky end 2</th>
                <th>Score</th>
                </tr>
            """
            for i in range(min(display_row,len(dna12_score))):
                print """
                    <tr>
                    <td>%.0f to %.0f</td>
                    <td>%s</td>
                    <td>%.0f to %.0f</td>
                    <td>%s</td>
                    <td>%.2f</td>
                    </tr>
                """ %(dna12_sorted[i,0], dna12_sorted[i,1], num_to_base(dna12_sorted[i,2:2+len(dna1_se[0,:])]),
                      dna12_sorted[i,2+len(dna1_se[0,:])], dna12_sorted[i,3+len(dna1_se[0,:])],
                      num_to_base(dna12_sorted[i,4+len(dna1_se[0,:]):4+len(dna1_se[0,:])+len(dna2_se[0,:])]),
                      dna12_sorted[i,-1])
            print '</table><br/>'
        
      
def input_check():
    form = cgi.FieldStorage()
    inputok = False
    error_code = 0
    output_code = 0
    se_score = 0
    display_row = int(form.getfirst("output_row"))
    if "dnaseq1" not in form:
        error_code = 0
    else:
        dna1_up, dna1_low, flag1 = bases_to_num(form.getfirst("dnaseq1"))
        if flag1 == 0:
            error_code = 1
        else:
            inputok = True
            dna_num = 1

            if "dnaseq2" in form:
                dna_num = 0
                dna2_up, dna2_low, flag2 = bases_to_num(form.getfirst("dnaseq2"))
                if flag2 == 0:
                    inputok = False
                    error_code = 2
                else:
                    inputok = True
                    dna_num = 2

    if inputok:
        if ((dna_num == 1) and ("se1_type" not in form)) or ((dna_num == 2) and (("se1_type" not in form) or ("se2_type" not in form))):
            inputok = False
            error_code = 3

    if inputok and ("delG_lig" in form) and (not is_float(form.getfirst("delG_lig"))):
        inputok = False
        error_code = 4
    elif ("delG_lig" in form) and (is_float(form.getfirst("delG_lig"))):
        delG_lig = float(form.getfirst("delG_lig"))
    else:
        delG_lig = -6.3
        
    if inputok:
        overhang1 = int(form.getfirst("overhang1"))
        if (len(dna1_up) < overhang1 + 2):
            inputok = False
            error_code = 5
        else: 
            dna1_se = numpy.zeros((len(dna1_up)-overhang1-1, overhang1))
            dna1_delG = numpy.zeros((len(dna1_up)-overhang1-1,6))
            dna1_score = numpy.zeros(len(dna1_up)-overhang1-1)
            for i in range(len(dna1_up)-overhang1-1):
                se1l = form.getfirst("dnaseq1")[i]    # keep as text
                se1_up = dna1_up[i+1:i+overhang1+1]
                se1_low = dna1_low[len(dna1_low)-overhang1-i-1:len(dna1_low)-i-1]
                se1r = form.getfirst("dnaseq1")[i+overhang1+1]    # keep as text
                delG1 = delG_cal(se1_up, se1_low, se1l, se1r, form.getfirst("se1_type"),
                             form.getfirst("stacking_opt"), delG_lig)
                delG1_mm = delG_mm_cal(se1_up, se1_low, se1l, se1r, form.getfirst("se1_type"),
                                form.getfirst("stacking_opt"), delG_lig)
                dna1_se[i] = se1_up
                dna1_delG[i,0] = delG1
                dna1_delG[i,1:6] = delG1_mm
                dna1_score[i] = score_cal(delG1, delG1_mm)

             # need to sort score

            if dna_num == 1:
                output_code = 0
                se_score = (dna1_se, dna1_score)
            elif dna_num == 2:
                overhang2 = int(form.getfirst("overhang2"))
                if (len(dna2_up) < overhang2 + 2):
                    inputok = False
                    error_code = 5
                else:
                    dna2_se = numpy.zeros((len(dna2_up)-overhang2-1, overhang2))
                    dna2_delG = numpy.zeros((len(dna2_up)-overhang2-1, 6))
                    dna2_score = numpy.zeros(len(dna2_up)-overhang2-1)
                    for i in range(len(dna2_up)-overhang2-1):
                        se2l = form.getfirst("dnaseq2")[i]
                        se2_up = dna2_up[i+1:i+overhang2+1]
                        se2_low = dna2_low[len(dna2_low)-overhang2-i-1:len(dna2_low)-i-1]
                        se2r = form.getfirst("dnaseq2")[i+overhang2+1]
                        delG2 = delG_cal(se2_up, se2_low, se2l, se2r, form.getfirst("se2_type"),
                                         form.getfirst("stacking_opt"), delG_lig)
                        delG2_mm = delG_mm_cal(se2_up, se2_low, se2l, se2r, form.getfirst("se2_type"),
                                               form.getfirst("stacking_opt"), delG_lig)
                        dna2_se[i] = se2_up
                        dna2_delG[i,0] = delG2
                        dna2_delG[i,1:6] = delG2_mm
                        dna2_score[i] = score_cal(delG2, delG2_mm)

                    if form.getfirst("se1_type") != form.getfirst("se2_type"):
                        output_code = 1
                        se_score = (dna1_se, dna1_score, dna2_se, dna2_score)
                    else:
                        output_code = 2
                        dna12_delG = numpy.zeros((len(dna1_se)*len(dna2_se),8))
                        dna12_score = numpy.zeros(len(dna1_se)*len(dna2_se))
                        for i in range(len(dna1_se)):
                            se1l = form.getfirst("dnaseq1")[i]    # keep as text
                            se1_up = dna1_up[i+1:i+overhang1+1]
                            se1_low = dna1_low[len(dna1_low)-overhang1-i-1:len(dna1_low)-i-1]
                            se1r = form.getfirst("dnaseq1")[i+overhang1+1]    # keep as text
                            for j in range(len(dna2_se)):
                                se2l = form.getfirst("dnaseq2")[j]
                                se2_up = dna2_up[j+1:j+overhang2+1]
                                se2_low = dna2_low[len(dna2_low)-overhang2-j-1:len(dna2_low)-j-1]
                                se2r = form.getfirst("dnaseq2")[j+overhang2+1]
                                delG1and2_mm = delG_2se_mm_cal(se1_up, se1_low, se2_up, se2_low,
                                                    se1l, se1r, se2l, se2r, form.getfirst("se1_type"),
                                                            form.getfirst("stacking_opt"), delG_lig)
                                dna12_delG[i*len(dna2_se) + j,:] = delG1and2_mm
                                dna12_score[i*len(dna2_se) + j] = score_cal_2se(dna1_delG[i,0], dna1_delG[i,1:6], dna2_delG[j,0], dna2_delG[j,1:6], delG1and2_mm)
                        se_score = (dna1_se, dna1_score, dna2_se, dna2_score, dna12_score)

                    
    return inputok, error_code, output_code, se_score, display_row

def output_html(inputok, error_code, output_code, se_score, display_row):
    print "Content-Type: text/html\n"
    print """
        <!DOCTYPE html>
        <html>
        <head>

        <title>NP-Sticky</title>
        <link rel="stylesheet" type="text/css" href="../seas/homestyle.css">

        </head>

        <body>
        <div id="content">

        <h1><a href="../seas/index.html" style="background-color:#C0C0C0; color:black;">NP-Sticky: Model-Guided Ligation Calculator</a></h1>

        <ul>
            <li><a href="../seas/about.html" >Overview</a></li>
            <li><a href="../seas/enzcompare.html">Type IIP vs. IIS Restriction Enzymes</a></li>
            <li><a href="../seas/schematic.html">Model Schematic</a></li>
            <li><a href="../seas/modelresults.html">Model Recommendations</a></li>
        </ul>
        <ul>
            <li><a href="../seas/stickyfinder.html" style="background-color:darkblue">Sticky End Finder</a></li>
            <li><a href="../seas/stickyend.html">Sticky End Free Energy Calculator</a></li>
        </ul>
        <ul>
            <li><a href="../seas/linligation.html">Linear Ligation Calculator</a></li>
            <li><a href="../seas/cirligation.html">Circular Ligation Calculator</a></li>
        </ul>

        <hr>

        <div id="white">

        <h2>Sticky End Finder - Results</h2>
                
        <p class="justify">
        """
    if not inputok:
        error_output(error_code)
    else:
        good_output(output_code, se_score, display_row)
    print """
        </p>

        </div>

        <hr>

        <p class="name">Last modified by Daphne Ng on Jan 4th, 2014</p>

        </div>
        </body>
        </html>
        """

           
inputok, error_code, output_code, se_score, display_row = input_check()
output_html(inputok, error_code, output_code, se_score, display_row)
