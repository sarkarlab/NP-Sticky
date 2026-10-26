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

def func_lin_2piece(x, *args):
    #x = A, B, AB
    K, conc = args
    y = [x[0]*x[1] - K*x[2],
         x[0] + x[2] - conc[0],
         x[1] + x[2] - conc[1]]
    return y

def solve_func(lib, *param):
    K, conc = param

    xini = [random.uniform(0, max(conc)) for i in range(3)]
    xans, dic, ier, mesg = fsolve(func_lin_2piece, xini, args=param, full_output=1)
    count = 0
    while ((ier != 1) or (min(xans) < 0)) and count < 50:
        count = count + 1
        xini = [random.uniform(0, max(conc)) for i in range(3)]
        xans, dic, ier, mesg = fsolve(func_lin_2piece, xini, args=param, full_output=1)
    incorp = (xans[2]/conc[lib-1])*100
    eff = (xans[2]/(conc[0]+conc[1]))*2*100

    return xans, incorp, eff, count

def compute_ll_recom(K, lib_conc, lib_opt):
    n = 30
    ratio = numpy.logspace(-3, 3, num=n, base=2.0)
    lib_range = numpy.linspace(lib_conc/10.0, lib_conc, num=n)    
    incorp, eff = numpy.zeros((n,n)), numpy.zeros((n,n))
               
    for j in range(n):
        for i in range(n):
            if lib_opt == 1:
                conc_param = [lib_range[j], ratio[i]*lib_range[j]]
            elif lib_opt == 2:
                conc_param = [ratio[i]*lib_range[j], lib_range[j]]
        
            out_conc, incorp[j,i], eff[j,i], count = solve_func(lib_opt, K, conc_param)
            if count >= 50:
                break

    fig_lin_incorp = plt.figure()
    plt.pcolormesh(ratio, lib_range, incorp)
    plt.title("Percent Incorporation")
    plt.xlabel("Molar ratio of [flanking DNAs]/[library]")
    plt.ylabel("Library concentration (nM)")
    plt.xscale("log", basex=2.0)
    plt.xticks((1.0/8, 1.0/4, 1.0/2, 1, 2, 4, 8), ("1/8", "1/4", "1/2", "1", "2", "4", "8"))
    plt.xlim((ratio[0], ratio[-1]))
    plt.ylim((lib_range[0], lib_range[-1]))
    plt.colorbar()
    outputobj1 = StringIO.StringIO()
    fig_lin_incorp.savefig(outputobj1, format="png")
    data_lin_incorp = outputobj1.getvalue().encode('base64').replace('\n', '')
    outputobj1.close()
    lincorp_image_tag = '<img src="data:image/png;base64,{0}" height="400" width="500">'.format(data_lin_incorp)
    

    fig_lin_eff = plt.figure()
    plt.pcolormesh(ratio, lib_range, eff)
    plt.title("Percent Efficiency")
    plt.xlabel("Molar ratio of [flanking DNAs]/[library]")
    plt.ylabel("Library concentration (nM)")
    plt.xscale("log", basex=2.0)
    plt.xticks((1.0/8, 1.0/4, 1.0/2, 1, 2, 4, 8), ("1/8", "1/4", "1/2", "1", "2", "4", "8"))
    plt.xlim((ratio[0], ratio[-1]))
    plt.ylim((lib_range[0], lib_range[-1]))
    plt.colorbar()
    outputobj2 = StringIO.StringIO()
    fig_lin_eff.savefig(outputobj2, format="png")
    data_lin_eff = outputobj2.getvalue().encode('base64').replace('\n', '')
    outputobj2.close()
    leff_image_tag = '<img src="data:image/png;base64,{0}" height="400" width="500">'.format(data_lin_eff)

    return lincorp_image_tag, leff_image_tag, count   
    
def error_output(flag):
    print "Error:"
    if flag == 0:
        print "Please indicate the DNA concentration type (molar concentration or mass concentration)."
    elif flag == 1:
        print "Please indicate the type of sticky end information."
    elif flag == 2:
        print "Please enter DNA concentration (if mass concentration is entered, please enter the number of base pairs)."
    elif flag == 3:
           print "Sticky end information missing."
    elif flag == 4:
        print "DNA concentration information invalid."
    elif flag == 5:
        print "Sticky end information invalid."
    elif flag == 6:
        print "Product distribution cannot be computed. Please check input values."

def good_output(lincorp_tag, leff_tag, lib_s):
    print "Percent incorporation is defined as [full-length product]/[library (DNA %s)]<sub>input</sub>.<br/>" %lib_s
    print lincorp_tag
    print "<br/><br/>"
    print "Percent efficiency is defined as 2 * [full-length product]/[all reactants]<sub>input</sub>.<br/>" ## need to change for 3
    print leff_tag
       
    
def import_check():
    form = cgi.FieldStorage()
    inputok = False
    error_code = 0
    lin_tag, leff_tag = "", ""
    
    if "conctype_lib_opt2" not in form:
        error_code = 0
    elif "se_opt2" not in form:
        error_code = 1
    else:
        if ((form.getfirst("conctype_lib_opt2") == "molar") and ("ic_lib_opt2" not in form)) or (
            (form.getfirst("conctype_lib_opt2") == "mass") and (("imc_lib_opt2" not in form) or ("bp_lib_opt2" not in form))):
            error_code = 2
        elif ((form.getfirst("se_opt2") == "delG") and ("delG12_opt2" not in form)) or (
              (form.getfirst("se_opt2") == "base") and ("se12_opt2" not in form)):
            error_code = 3
        else:
            inputok = True

    if inputok:
        inputok = False
        conctype_lib = form.getfirst("conctype_lib_opt2")
        se_type = form.getfirst("se_opt2")

        if ((conctype_lib == "molar") and not (is_float(form.getfirst("ic_lib_opt2")))) or (
            (conctype_lib == "mass") and not (is_float(form.getfirst("imc_lib_opt2")) or (is_integer(form.getfirst("bp_lib_opt2"))))):
            error_code = 4
        elif ((se_type == "delG") and not (is_float(form.getfirst("delG12_opt2")))):
            error_code = 5
        elif (se_type == "base"):
            se12_up, se12_low, flag = bases_to_num(form.getfirst("se12_opt2"))
            if (flag == 0):
                error_code = 5
            else:
                inputok = True
        else:
            inputok = True

    if inputok:
        if conctype_lib == "molar":
            ic_lib = float(form.getfirst("ic_lib_opt2"))
        elif conctype_lib == "mass":
            ic_lib = mass_to_molar_stock(form.getfirst("imc_lib_opt2"), form.getfirst("bp_lib_opt2"))

        if se_type == "delG":
            K = delG_to_K(float(form.getfirst("delG12_opt2")))
        elif se_type == "base":
            delG_lig = -6.3    #kcal/mol
            delG = delG_cal(se12_up, se12_low, form.getfirst("se12_l_opt2"),
                     form.getfirst("se12_r_opt2"), form.getfirst("se12_type_opt2"), delG_lig)
            K = delG_to_K(delG)                        

        libopt = int(form.getfirst("lib_opt2"))
        lin_tag, leff_tag, count = compute_ll_recom(K, ic_lib, libopt)
        if count >= 50:
            error_code = 6
            inputok = False
        
    return inputok, error_code, lin_tag, leff_tag, form.getfirst("lib_opt2")

def output_html(inputok, error_code, lin_tag, leff_tag, lib_s):
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
            <li><a href="../seas/linligation.html" style="background-color:darkblue">Linear Ligation Calculator</a></li>
            <li><a href="../seas/cirligation.html">Circular Ligation Calculator</a></li>
        </ul>

        <hr>

        <div id="white">

        <h2>Linear Ligation Calculator - Results</h2>
        <h3>Ligation Recommendation</h3>
        
        <p class="justify">
        """
    if not inputok:
        error_output(error_code)
    else:
        good_output(lin_tag, leff_tag, lib_s)
    print """
        </p>

        </div>

        <hr>

        <p class="name">Last modified by Daphne Ng on Jan 4th, 2014</p>

        </div>
        </body>
        </html>
        """
    

inputok, error_code, lin_tag, leff_tag, lib_s= import_check()
output_html(inputok, error_code, lin_tag, leff_tag, lib_s)
