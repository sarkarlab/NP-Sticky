#!/usr/bin/python

#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License.

import cgi

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
        print "Information of sticky end 1 missing."
    elif flag == 1:
        print "Please only enter A/C/T/G for sticky end 1."
    elif flag == 2:
        print "Information of sticky end 2 missing."
    elif flag == 3:
        print "Please only enter A/C/T/G for sticky end 2."
    elif flag == 4:
        print "Invalid custom free energy for the repair of phosphodiester bonds."


def good_output(opt, delGs):
    if opt == 0:
        delG1, delG1_mm = delGs
    elif opt == 1:
        delG1, delG1_mm, delG2, delG2_mm = delGs
    elif opt == 2:
        delG1, delG1_mm, delG2, delG2_mm, delG1and2_mm = delGs
    uu, ll, uu_gap, ll_gap, ul_gap = delG1_mm
    print """
        <b>Sticky end 1 (SE1):</b><br/>
        The free energy of sticky end 1 is <b>%.2f kcal/mol</b>.<br/>
        <table border="1">
        <tr>
        <th>Sticky end mismatch</th>
        <th>Free energy (kcal/mol)</th>
        </tr>
        <tr>
        <td>SE1 upper strand - SE1 upper strand</td>
        <td>%.2f
        <tr/>
        <tr>
        <td>SE1 lower strand - SE1 lower strand</td>
        <td>%.2f
        </tr>
        <tr>
        <td>SE1 upper strand - SE1 upper strand (with 1 nucleotide gap)</td>
        <td>%.2f
        </tr>
        <tr>
        <td>SE1 lower strand - SE1 lower strand (with 1 nucleotide gap)</td>
        <td>%.2f
        </tr>
        <tr>
        <td>SE1 upper strand - SE1 lower strand (with 1 nucleotide gap)</td>
        <td>%.2f
        </tr>
        </table>
        <br/><br/>
    """ % (delG1, uu, ll, uu_gap, ll_gap, ul_gap)

    if opt >= 1:
        uu, ll, uu_gap, ll_gap, ul_gap = delG2_mm
        print """
            <b>Sticky end 2 (SE2):</b><br/>
            The free energy of sticky end 2 is <b>%.2f kcal/mol</b>.<br/><br/>
            <table border="1">
            <tr>
            <th>Sticky end mismatch</th>
            <th>Free energy (kcal/mol)</th>
            </tr>
            <tr>
            <td>SE2 upper strand - SE2 upper strand</td>
            <td>%.2f
            <tr/>
            <tr>
            <td>SE2 lower strand - SE2 lower strand</td>
            <td>%.2f
            </tr>
            <tr>
            <td>SE2 upper strand - SE2 upper strand (with 1 nucleotide gap)</td>
            <td>%.2f
            </tr>
            <tr>
            <td>SE2 lower strand - SE2 lower strand (with 1 nucleotide gap)</td>
            <td>%.2f
            </tr>
            <tr>
            <td>SE2 upper strand - SE2 lower strand (with 1 nucleotide gap)</td>
            <td>%.2f
            </tr>
            </table>
            <br/><br/>
        """ % (delG2, uu, ll, uu_gap, ll_gap, ul_gap)
		
        print "<b>Mismatch of sticky end 1 (SE1) with sticky end 2 (SE2):</b><br/>"

        if opt == 1:
            print "Different types of overhangs, no mismatch!"
        elif opt == 2:
            u1u2, l1l2, u1l2, l1u2, u1u2_gap, l1l2_gap, u1l2_gap, l1u2_gap = delG1and2_mm
            print """
                <table border="1">
                <tr>
                <th>Sticky end mismatch</th>
                <th>Free energy (kcal/mol)</th>
                </tr>
                <tr>
                <td>SE1 upper strand - SE2 upper strand</td>
                <td>%.2f
                <tr/>
                <tr>
                <td>SE1 lower strand - SE2 lower strand</td>
                <td>%.2f
                </tr>
                <tr>
                <td>SE1 upper strand - SE2 lower strand</td>
                <td>%.2f
                <tr/>
                <tr>
                <td>SE1 lower strand - SE2 upper strand</td>
                <td>%.2f
                </tr>
                <tr>
                <td>SE1 upper strand - SE2 upper strand (with 1 nucleotide gap)</td>
                <td>%.2f
                </tr>
                <tr>
                <td>SE1 lower strand - SE2 lower strand (with 1 nucleotide gap)</td>
                <td>%.2f
                </tr>
                <tr>
                <td>SE1 upper strand - SE2 lower strand (with 1 nucleotide gap)</td>
                <td>%.2f
                </tr>
                <tr>
                <td>SE1 lower strand - SE2 upper strand (with 1 nucleotide gap)</td>
                <td>%.2f
                </tr>
                </table>
            """ % (u1u2, l1l2, u1l2, l1u2, u1u2_gap, l1l2_gap, u1l2_gap, l1u2_gap)
            

        
def input_check():
    form = cgi.FieldStorage()
    inputok = False
    error_code = 0
    output_code = 0
    delGs = 0
    if "se1" not in form or "se1_type" not in form:
        error_code = 0
    else:
        sticky1_up, sticky1_low, flag1 = bases_to_num(form.getfirst("se1"))
        if flag1 == 0:
            error_code = 1
        else:
            inputok = True
            se_num = 1

            if "se2" in form:
                se_num = 0
                if "se2_type" not in form:
                    inputok = False
                    error_code = 2
                else:
                    sticky2_up, sticky2_low, flag2 = bases_to_num(form.getfirst("se2"))
                    if flag2 == 0:
                        inputok = False
                        error_code = 3
                    else:
                        inputok = True
                        se_num = 2

    if inputok and ("delG_lig" in form) and (not is_float(form.getfirst("delG_lig"))):
        inputok = False
        error_code = 4
    elif "delG_lig" in form:
        delG_lig = float(form.getfirst("delG_lig"))
    else:
        delG_lig = -6.3
        
    if inputok:
        delG1 = delG_cal(sticky1_up, sticky1_low, form.getfirst("se1l"),
                     form.getfirst("se1r"), form.getfirst("se1_type"),
                         form.getfirst("stacking_opt"), delG_lig)
        delG1_mm = delG_mm_cal(sticky1_up, sticky1_low, form.getfirst("se1l"),
                        form.getfirst("se1r"), form.getfirst("se1_type"),
                            form.getfirst("stacking_opt"), delG_lig)
        if se_num == 1:
            output_code = 0
            delGs= (delG1, delG1_mm)
        elif se_num == 2:
            delG2 = delG_cal(sticky2_up, sticky2_low, form.getfirst("se2l"),
                         form.getfirst("se2r"), form.getfirst("se2_type"),
                             form.getfirst("stacking_opt"), delG_lig)
            delG2_mm = delG_mm_cal(sticky2_up, sticky2_low, form.getfirst("se2l"),
                            form.getfirst("se2r"), form.getfirst("se2_type"),
                                form.getfirst("stacking_opt"), delG_lig)
            if form.getfirst("se1_type") != form.getfirst("se2_type"):
                output_code = 1
                delGs = (delG1, delG1_mm, delG2, delG2_mm)
            else:    #also notice mismatch score for different lengths of delG is problemmatics
                delG1and2_mm = delG_2se_mm_cal(sticky1_up, sticky1_low, sticky2_up, sticky2_low,
                                    form.getfirst("se1l"), form.getfirst("se1r"), form.getfirst(
                                        "se2l"), form.getfirst("se2r"), form.getfirst("se1_type"),
                                            form.getfirst("stacking_opt"), delG_lig)
                output_code = 2
                delGs = (delG1, delG1_mm, delG2, delG2_mm, delG1and2_mm)

    return inputok, error_code, output_code, delGs


def output_html(inputok, error_code, output_code, delGs):
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
            <li><a href="../seas/stickyend.html" style="background-color:darkblue">Sticky End Free Energy Calculator</a></li>
        </ul>
        <ul>
            <li><a href="../seas/linligation.html">Linear Ligation Calculator</a></li>
            <li><a href="../seas/cirligation.html">Circular Ligation Calculator</a></li>
        </ul>

        <hr>

        <div id="white">

        <h2>Sticky End Free Energy Calculator - Results</h2>
        <h3>By DNA Thermodynamics</h3>
        
        <p class="justify">
        """
    if not inputok:
        error_output(error_code)
    else:
        good_output(output_code, delGs)
    print """
        </p>
        <br/>
        <p class="justify"><i><b>Terminology:</i></b></p>
        <p class="center"><img src="../seas/figures/web_SEcalc.png" type="image/png" width="309" height="115"></p>
        
        </div>

        <hr>

        <p class="name">Last modified by Daphne Ng on Jan 14th, 2014</p>

        </div>
        </body>
        </html>
        """

           
inputok, error_code, output_code, delGs = input_check()
output_html(inputok, error_code, output_code, delGs)
