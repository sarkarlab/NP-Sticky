#!/usr/bin/python

#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License.

import cgi

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy
import random
from scipy.optimize import fsolve

import StringIO

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def is_integer(s):
    try:
        int(s)
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

def mass_to_molar_stock(mass_conc_s, bp_s):
    mass_conc = float(mass_conc_s)
    bp = float(bp_s)
    MW = bp*607.4 + 157.9
    molar_conc = (mass_conc*0.001)/MW*1000000000
    return molar_conc

def delG_to_K(delG):
    T = 25 + 273.15
    R = 1.9858775
    K = numpy.exp((delG*1e3)/(R*T))*1e9    #nM
    return K

def Kcir_cal(K, bp):
    MW = bp*607.4 + 157.9    #g/mol
    j = 3.6*((30.8e6/MW)**(1.5))    #*1e11 ends/ml-original eq 3.6e11*((30.8e6/MW)**(1.5))
    Aeff = j/6.022    #nM - original eq (j/6.022e23)*1e3*1e9
    Kcir = K/Aeff    # no units
    return Kcir

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
        
def stack(left1, left2, right1, right2):
    STACKSCORE = [[-0.8, -0.5, -0.3, -0.9],
                      [-0.5, -0.4, -0.2, -0.5],
                      [-0.3, -0.2, -0.6, -0.6],
                      [-0.9, -0.5, -0.6, -1.1]]
    return STACKSCORE[left1][left2] + STACKSCORE[right1][right2]

def delG_cal(se_up, se_low, se_left, se_right, se_type, delG_lig):
    se_left_up, se_left_low, flag = bases_to_num(se_left)
    se_right_up, se_right_low, flag = bases_to_num(se_right)
    if se_type == "5":
        delG_se = h_bond(se_up, se_low) + stack(se_left_up[0], se_up[0], se_right_low[0], se_low[0]) + delG_lig
    elif se_type == "3":
        delG_se = h_bond(se_up, se_low) + stack(se_low[-1], se_left_low[0], se_up[-1], se_right_up[0]) + delG_lig

    return delG_se

def func_cir_VI(x, *args):
    #x = V, I, VI1, VI2, VIc1, VIc2, IVI1, IVI2, VIV1, VIV2,
    #    VIVI1a, VIVI2a, VIVI1b, VIVI2b, VIVIc1, VIVIc2

    Ki, Kcir_VI, Kcir_VIVI, conc = args    #conc = V, I
    y = [x[0]*x[1] - Ki[0]*x[2],
         x[0]*x[1] - Ki[1]*x[3],
         x[2] - Kcir_VI[1]*x[5],
         x[3] - Kcir_VI[0]*x[4],
         x[2]*x[1] - Ki[1]*x[7],
         x[3]*x[1] - Ki[0]*x[6],
         x[2]*x[0] - Ki[1]*x[9],
         x[3]*x[0] - Ki[0]*x[8],
         (x[6]+x[7])*x[0] - Ki[0]*x[10],
         (x[6]+x[7])*x[0] - Ki[1]*x[11],
         (x[8]+x[9])*x[1] - Ki[0]*x[12],
         (x[8]+x[9])*x[1] - Ki[1]*x[13],
         (x[10]+x[12]) - Kcir_VIVI[1]*x[15],
         (x[11]+x[13]) - Kcir_VIVI[0]*x[14],
         x[0]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+2*(x[8]+x[9]
            +x[10]+x[11]+x[12]+x[13]+x[14]+x[15]) - conc[0],
         x[1]+x[2]+x[3]+x[4]+x[5]+x[8]+x[9]+2*(x[6]+x[7]
            +x[10]+x[11]+x[12]+x[13]+x[14]+x[15]) - conc[1]]
    return y

def cir_solve_func(*param):
    Ki, Kcir_VI, Kcir_VIVI, conc = param
    xini = [random.uniform(0,max(conc)) for i in range(16)]
    xans, dic, ier, mesg = fsolve(func_cir_VI, xini, args=param, full_output=1)
    count = 0
    while ((ier != 1) or (min(xans) < 0)) and count < 50:
        count = count + 1
        xini = [random.uniform(0,max(conc)) for i in range(16)]
        xans, dic, ier, mesg = fsolve(func_cir_VI, xini, args=param, full_output=1)
    incorp = (xans[4]+xans[5])/conc[1]*100
    eff = (xans[4]+xans[5])/(conc[0]+conc[1])*2*100
    accur = (xans[4]+xans[5])/(xans[4]+xans[5]+xans[14]+xans[15])*100

    return xans, incorp, eff, accur, count
'''
def format_input(ic_lib_s, delG_s, bp_s):
    T = 25 + 273.15
    R = 1.9858775

    delG = numpy.zeros(2)
    for i in range(2):
        delG[i] = float(delG_s[i])
    ic_lib = float(ic_lib_s)
    bp = numpy.array([int(bp_s), 2*int(bp_s)])

    Ki = numpy.exp(delG*1e3/(R*T))*1e9    #nM
    MW = bp*607.4 + 157.9    #g/mol
    j = 3.6*((30.8e6/MW)**(1.5))    #*1e11 ends/ml-original eq 3.6e11*((30.8e6/MW)**(1.5))
    Aeff = j/6.022    #nM - original eq (j/6.022e23)*1e3*1e9
    Kcir_VI = Ki/Aeff[0]    # no units
    Kcir_VIVI = Ki/Aeff[1]    # no units

    return Ki, Kcir_VI, Kcir_VIVI, ic_lib
'''
def compute_cl_recom(Ki, Kcir_VI, Kcir_VIVI, lib_conc):
    n = 30
    VIratio_range = numpy.logspace(-3, 3, num=n, base=2.0)
    lib_conc_range = numpy.linspace(lib_conc/10.0, lib_conc, num=n)
    incorp, eff, accur = numpy.zeros((n,n)), numpy.zeros((n,n)), numpy.zeros((n,n))

    for j in range(n):
        for i in range(n):
            out_conc, incorp[j,i], eff[j,i], accur[j,i], count = cir_solve_func(Ki,
                Kcir_VI, Kcir_VIVI, [VIratio_range[i]*lib_conc_range[j], lib_conc_range[j]])
            if count >= 50:
                break

    fig_cir_incorp = plt.figure()
    plt.pcolormesh(VIratio_range, lib_conc_range, incorp)
    plt.title("Percent Incorporation")
    plt.xlabel("Molar ratio of [vector]/[insert library]")
    plt.ylabel("Insert concentration (nM)")
    plt.xscale("log", basex=2.0)
    plt.xticks((1.0/8, 1.0/4, 1.0/2, 1, 2, 4, 8), ("1/8", "1/4", "1/2", "1", "2", "4", "8"))
    plt.xlim((VIratio_range[0], VIratio_range[-1]))
    plt.ylim((lib_conc_range[0], lib_conc_range[-1]))
    plt.colorbar()
    outputobj1 = StringIO.StringIO()
    fig_cir_incorp.savefig(outputobj1, format="png")
    data_cir_incorp = outputobj1.getvalue().encode('base64').replace('\n', '')
    outputobj1.close()
    cincorp_image_tag = '<img src="data:image/png;base64,{0}" height="400" width="500">'.format(data_cir_incorp)

    fig_cir_eff = plt.figure()
    plt.pcolormesh(VIratio_range, lib_conc_range, eff)
    plt.title("Percent Efficiency")
    plt.xlabel("Molar ratio of [vector]/[insert library]")
    plt.ylabel("Insert concentration (nM)")
    plt.xscale("log", basex=2.0)
    plt.xticks((1.0/8, 1.0/4, 1.0/2, 1, 2, 4, 8), ("1/8", "1/4", "1/2", "1", "2", "4", "8"))
    plt.xlim((VIratio_range[0], VIratio_range[-1]))
    plt.ylim((lib_conc_range[0], lib_conc_range[-1]))
    plt.colorbar()
    outputobj2 = StringIO.StringIO()
    fig_cir_eff.savefig(outputobj2, format="png")
    data_cir_eff = outputobj2.getvalue().encode('base64').replace('\n', '')
    outputobj2.close()
    ceff_image_tag = '<img src="data:image/png;base64,{0}" height="400" width="500">'.format(data_cir_eff)
    
    fig_cir_acc = plt.figure()
    plt.pcolormesh(VIratio_range, lib_conc_range, accur)
    plt.title("Percent Accuracy")
    plt.xlabel("Molar ratio of [vector]/[insert library]")
    plt.ylabel("Insert concentration (nM)")
    plt.xscale("log", basex=2.0)
    plt.xticks((1.0/8, 1.0/4, 1.0/2, 1, 2, 4, 8), ("1/8", "1/4", "1/2", "1", "2", "4", "8"))
    plt.xlim((VIratio_range[0], VIratio_range[-1]))
    plt.ylim((lib_conc_range[0], lib_conc_range[-1]))     
    plt.colorbar()
    
    outputobj3 = StringIO.StringIO()
    fig_cir_acc.savefig(outputobj3, format="png")
    data_cir_acc = outputobj3.getvalue().encode('base64').replace('\n', '')
    outputobj3.close()
    cacc_image_tag = '<img src="data:image/png;base64,{0}" height="400" width="500">'.format(data_cir_acc)

    return cincorp_image_tag, ceff_image_tag, cacc_image_tag, count

def error_output(flag):
    print "Error:"
    if flag == 0:
        print "Plesae enter the number of base pairs for both insert and vector."
    elif flag == 1:
        print "Please indicate the DNA concentration type (molar concentration or mass concentration)."
    elif flag == 2:
        print "Please enter DNA concentration."
    elif flag == 3:
        print "Please indicate the type of sticky end information."
    elif flag == 4:
        print "Sticky end information missing."
    elif flag == 5:
        print "Base pair information invalid."
    elif flag == 6:
        print "DNA concentration information invalid."
    elif flag == 7:
        print "Sticky end information invalid."
    elif flag == 8:
        print "Product distribution cannot be computed. Please check input values."

def good_output(cin_tag, ceff_tag, cacc_tag):
    print "Percent incorporation is defined as [circular vector-insert product]/[insert concentration]<sub>input</sub>.<br/>"
    print cin_tag
    print "<br/><br/>"
    print "Percent efficiency is defined as [circular vector-insert product]/[vector and insert]<sub>input</sub>.<br/>"
    print ceff_tag
    print "<br/><br/>"
    print "Percent accuracy is defined as [circular vector-insert product]/[all circular products].<br/>"
    print cacc_tag


def input_check():
    form = cgi.FieldStorage()
    inputok = False
    error_code = 0
    cin_tag, ceff_tag, cacc_tag = "", "", ""
    if ("bp_i" not in form) or ("bp_v" not in form):
        error_code = 0
    elif ("conctype_i" not in form):
        error_code = 1
    else:
        if ((form.getfirst("conctype_i") == "molar") and ("in_conc_i" not in form)) or (
            (form.getfirst("conctype_i") == "mass") and ("in_mconc_i" not in form)):
            error_code = 2
        else:
            if "se_opt" not in form:
                error_code = 3
            else:
                if ((form.getfirst("se_opt") == "delG") and (("delG12" not in form) or ("delG23" not in form))) or (
                        (form.getfirst("se_opt") == "base") and (("se12" not in form) or ("se23" not in form))):
                        error_code = 4
                else:
                    inputok = True

    if inputok:
        inputok = False
        conctype_i = form.getfirst("conctype_i")
        se_opt = form.getfirst("se_opt")
        if not (is_integer(form.getfirst("bp_i"))) or not (is_integer(form.getfirst("bp_v"))):
            error_code = 5
        elif ((conctype_i == "molar") and not (is_float(form.getfirst("in_conc_i")))) or (
            (conctype_i == "mass") and (not (is_float(form.getfirst("in_mconc_i"))))):
            error_code = 6
        else:
            if (se_opt == "delG") and (not (is_float(form.getfirst("delG12"))) or not (is_float(form.getfirst("delG23")))):
                error_code = 7
            elif se_opt == "base":
                se12_up, se12_low, flag1 = bases_to_num(form.getfirst("se12"))
                se23_up, se23_low, flag2 = bases_to_num(form.getfirst("se23"))
                if (flag1 == 0) or (flag2 == 0):
                    error_code = 7
                else:
                    inputok = True
            else:
                inputok = True

    if inputok:
        if conctype_i == "molar":
            in_conc_i = float(form.getfirst("in_conc_i"))
        elif conctype_i == "mass":
            in_conc_i = mass_to_molar_stock(form.getfirst("in_mconc_i"), form.getfirst("bp_i"))

        if se_opt == "delG":
            K12 = delG_to_K(float(form.getfirst("delG12")))
            K23 = delG_to_K(float(form.getfirst("delG23")))
        elif se_opt == "base":
            delG_lig = -6.3
            delG12 = delG_cal(se12_up, se12_low, form.getfirst("se12_l"),
                     form.getfirst("se12_r"), form.getfirst("se12_type"), delG_lig)
            K12 = delG_to_K(delG12)
            delG23 = delG_cal(se23_up, se23_low, form.getfirst("se23_l"),
                     form.getfirst("se23_r"), form.getfirst("se23_type"), delG_lig)
            K23 = delG_to_K(delG23)

        K12cir_VI = Kcir_cal(K12, int(form.getfirst("bp_i")) + int(form.getfirst("bp_v")))
        K23cir_VI = Kcir_cal(K23, int(form.getfirst("bp_i")) + int(form.getfirst("bp_v")))
        K12cir_VIVI = Kcir_cal(K12, 2*(int(form.getfirst("bp_i")) + int(form.getfirst("bp_v"))))
        K23cir_VIVI = Kcir_cal(K23, 2*(int(form.getfirst("bp_i")) + int(form.getfirst("bp_v"))))
        cin_tag, ceff_tag, cacc_tag, count = compute_cl_recom((K12, K23), (K12cir_VI, K23cir_VI), (K12cir_VIVI, K23cir_VIVI), in_conc_i)
        if count >= 50:
            error_code = 8
            inputok = False
            
        '''
    if ("in_conc_i" not in form) or ("bp_iv" not in form) or ("delG_1" not
        in form) or ("delG_2" not in form):
        error_code = 0
    elif not (is_float(form.getfirst("in_conc_i")) and is_integer(
        form.getfirst("bp_iv")) and is_float(form.getfirst("delG_1"))
            and is_float(form.getfirst("delG_2"))):
        error_code = 1
    else:
        inputok = True
        Ki, Kcir_VI, Kcir_VIVI, ic_lib = format_input (form.getfirst(
            "in_conc_i"), [form.getfirst("delG_1"), form.getfirst("delG_2")],
                form.getfirst("bp_iv"))
        cin_tag, ceff_tag, cacc_tag, count = compute_cl_recom(Ki, Kcir_VI, Kcir_VIVI, ic_lib)
        if count >= 50:
            error_code = 2
            inputok = False
        '''
    return inputok, error_code, cin_tag, ceff_tag, cacc_tag

def output_html(inputok, error_code, cin_tag, ceff_tag, cacc_tag):
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
            <li><a href="../seas/stickyfinder.html">Sticky End Finder</a></li>
            <li><a href="../seas/stickyend.html">Sticky End Free Energy Calculator</a></li>
        </ul>
        <ul>
            <li><a href="../seas/linligation.html">Linear Ligation Calculator</a></li>
            <li><a href="../seas/cirligation.html" style="background-color:darkblue">Circular Ligation Calculator</a></li>
        </ul>

        <hr>

        <div id="white">

        <h2>Circular Ligation Calculator - Results</h2>
        <h3>Ligation Recommendation</h3>
        
        <p class="justify">
        """
    if not inputok:
        error_output(error_code)
    else:
        good_output(cin_tag, ceff_tag, cacc_tag)
    print """
        </p>

        </div>

        <hr>

        <p class="name">Last modified by Daphne Ng on Jan 4th, 2014</p>

        </div>
        </body>
        </html>
        """

inputok, error_code, cin_tag, ceff_tag, cacc_tag = input_check()
output_html(inputok, error_code, cin_tag, ceff_tag, cacc_tag)
