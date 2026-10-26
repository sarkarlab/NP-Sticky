#!/usr/bin/python

#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License.

import cgi

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy
import random
from scipy.optimize import fsolve

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

def mass_to_molar(mass_conc_s, bp_s, vol_s, total_vol_s):
    mass_conc = float(mass_conc_s)
    bp = float(bp_s)
    vol = float(vol_s)
    total_vol = float(total_vol_s)
    MW = bp*607.4 + 157.9
    molar_conc = ((mass_conc*vol)/MW)/(total_vol*0.000001)
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

def func_lin_nomm(x, *args):
    #x = A, B, AB
    K, conc = args
    y = [x[0]*x[1] - K*x[2],
         x[0] + x[2] - conc[0],
         x[1] + x[2] - conc[1]]
    return y

def solve_func(lib_s, *param):
    K, conc = param
    lib = int(lib_s)

    xini = [random.uniform(0, max(conc)) for i in range(3)]
    xans, dic, ier, mesg = fsolve(func_lin_nomm, xini, args=param, full_output=1)
    count = 0
    while ((ier != 1) or (min(xans) < 0)) and count < 50:
        count = count + 1
        xini = [random.uniform(0, max(conc)) for i in range(3)]
        xans, dic, ier, mesg = fsolve(func_lin_nomm, xini, args=param, full_output=1)
    incorp = (xans[2]/conc[lib-1])*100
    eff = (xans[2]/(conc[0]+conc[1]))*2*100
    return xans, incorp, eff, count



def error_output(flag):
    print "Error:"
    if flag == 0:
        print "Please select a valuable DNA fragment."
    elif flag == 1:
        print "Please indicate the DNA concentration type (molar concentration or mass concentration)."
    elif flag == 2:
        print "Please enter DNA concentration (if mass concentration is entered, please enter the number of base pairs, volume of DNA and total reaction volume)."
    elif flag == 3:
        print "Please indicate the type of sticky end information."
    elif flag == 4:
        print "Sticky end information missing."
    elif flag == 5:
        print "DNA concentration information invalid."
    elif flag == 6:
        print "Sticky end information invalid."
    elif flag == 7:
        print "Product distribution cannot be computed. Please check input values."


def good_output(out_conc, incorp, eff, lib_s):
    print """
        Incorporation is <b>%.1f percent</b> (defined as [product]/[DNA %s]<sub>input</sub>).<br/>
        Efficiency is <b>%.1f percent</b> (defined as 2 * [product]/[all reactants]<sub>input</sub>).<br/><br/>
        """ % (incorp, lib_s, eff)
    print """
        <table border="2">
            <tr>
                <th>DNA species</th>
                <th>DNA concentration (nM)</th>
            </tr>
            <tr>
                <td>DNA 1</td>
                <td>%.2f</td>
            </tr>
            <tr>
                <td>DNA 2</td>
                <td>%.2f</td>
            </tr>
            <tr>
                <td><b>Product: DNA 1 + 2</b></td>
                <td><b>%.2f</b></td>
            </tr>
        </table>
        """ % (out_conc[0], out_conc[1], out_conc[2])
        
    
def import_check():
    form = cgi.FieldStorage()
    inputok = False
    error_code, out_conc, incorp, eff = 0, 0, 0, 0
    if "lib_opt2" not in form:
        error_code = 0
    else:
        if ("conctype_dna1_opt2" not in form) or ("conctype_dna2_opt2" not in form):
            error_code = 1
        else:
            if ((form.getfirst("conctype_dna1_opt2") == "molar") and ("ic_dna1_opt2" not in form)) or (
                (form.getfirst("conctype_dna1_opt2") == "mass") and (("imc_dna1_opt2" not in form) or ("bp_dna1_opt2" not in form) or
                                                                     ("vol_dna1_opt2" not in form) or ("total_vol_opt2" not in form))):
                error_code = 2
            elif ((form.getfirst("conctype_dna2_opt2") == "molar") and ("ic_dna2_opt2" not in form)) or (
                  (form.getfirst("conctype_dna2_opt2") == "mass") and (("imc_dna2_opt2" not in form) or ("bp_dna2_opt2" not in form) or
                                                                     ("vol_dna2_opt2" not in form) or ("total_vol_opt2" not in form))):
                error_code = 2 
            else:
                if "se_opt2" not in form:
                    error_code = 3
                else:
                    if ((form.getfirst("se_opt2") == "delG") and ("delG12_opt2" not in form)) or (
                        (form.getfirst("se_opt2") == "base") and ("se12_opt2" not in form)):
                        error_code = 4
                    else:
                        inputok = True

    if inputok:
        inputok = False
        conctype_dna1 = form.getfirst("conctype_dna1_opt2")
        conctype_dna2 = form.getfirst("conctype_dna2_opt2")
        se_type = form.getfirst("se_opt2") 
        if ((conctype_dna1 == "molar") and not (is_float(form.getfirst("ic_dna1_opt2")))) or (
            (conctype_dna1 == "mass") and (not (is_float(form.getfirst("imc_dna1_opt2"))) or not (is_integer(form.getfirst("bp_dna1_opt2"))) or
                                           not (is_float(form.getfirst("vol_dna1_opt2"))) or not (is_float(form.getfirst("total_vol_opt2"))))):
            error_code = 5
        elif ((conctype_dna2 == "molar") and not (is_float(form.getfirst("ic_dna2_opt2")))) or (
              (conctype_dna2 == "mass") and (not (is_float(form.getfirst("imc_dna2_opt2"))) or not (is_integer(form.getfirst("bp_dna2_opt2"))) or
                                           not (is_float(form.getfirst("vol_dna2_opt2"))) or not (is_float(form.getfirst("total_vol_opt2"))))):
            error_code = 5
        else:
            if ((se_type == "delG") and not (is_float(form.getfirst("delG12_opt2")))):
                error_code = 6
            elif se_type == "base":
                se12_up, se12_low, flag = bases_to_num(form.getfirst("se12_opt2"))
                if (flag == 0):
                    error_code = 6
                else:
                    inputok = True
            else:
                inputok = True

    if inputok:
        if conctype_dna1 == "molar":
            ic_dna1 = float(form.getfirst("ic_dna1_opt2"))
        elif conctype_dna1 == "mass":
            ic_dna1 = mass_to_molar(form.getfirst("imc_dna1_opt2"), form.getfirst("bp_dna1_opt2"), form.getfirst("vol_dna1_opt2"), form.getfirst("total_vol_opt2"))
        if conctype_dna2 == "molar":
            ic_dna2 = float(form.getfirst("ic_dna2_opt2"))
        elif conctype_dna2 == "mass":
            ic_dna2 = mass_to_molar(form.getfirst("imc_dna2_opt2"), form.getfirst("bp_dna2_opt2"), form.getfirst("vol_dna2_opt2"), form.getfirst("total_vol_opt2"))

        if se_type == "delG":
            K = delG_to_K(float(form.getfirst("delG12_opt2")))
        elif se_type == "base":
            delG_lig = -6.3    #kcal/mol
            delG = delG_cal(se12_up, se12_low, form.getfirst("se12_l_opt2"),
                     form.getfirst("se12_r_opt2"), form.getfirst("se12_type_opt2"), delG_lig)
            K = delG_to_K(delG)

        out_conc, incorp, eff, count = solve_func(form.getfirst("lib_opt2"), K, (ic_dna1, ic_dna2))
        if count >= 50:
            error_code = 7
            inputok = False
            
    return inputok, error_code, out_conc, incorp, eff, form.getfirst("lib_opt2")

def output_html(inputok, error_code, out_conc, incorp, eff, lib_s):
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
        <h3>Product Distribution Prediction</h3>
        
        <p class="justify">
        """
    if not inputok:
        error_output(error_code)
    else:
        good_output(out_conc, incorp, eff, lib_s)
    print """
        </p>

        </div>

        <hr>

        <p class="name">Last modified by Daphne Ng on Jan 4th, 2014</p>

        </div>
        </body>
        </html>
        """

inputok, error_code, out_conc, incorp, eff, lib_s = import_check()
output_html(inputok, error_code, out_conc, incorp, eff, lib_s)
