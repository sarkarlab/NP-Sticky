#!/usr/bin/python

#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License.

import cgi

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy

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

def fitting_delG(trial, ic, omr, bp, g_low, g_high):
    T = 25 + 273.15
    R = 1.9858775

    g_space = numpy.arange(g_low, g_high + 0.01, 0.01)
    leastsq = numpy.zeros(trial)
    residue = numpy.zeros(len(g_space))
    j = 0

    for g in g_space:
        
        for i in range(trial):
            K = numpy.exp((g*1e3)/(R*T))*1e9    #nM
            AB = 0.5*((ic[i,0] + ic[i,1]+K) - numpy.sqrt((ic[i,0]
                       + ic[i,1]+K)**2 - 4*ic[i,0]*ic[i,1]))
            A = ic[i,0] - AB
            B = ic[i,1] - AB
            mol_ratio = numpy.array([A, B, AB])/AB
            bp_ratio = numpy.append(bp, sum(bp))/(sum(bp)*1.0)
            mass_ratio = mol_ratio * bp_ratio
            leastsq[i] = sum((mass_ratio - omr[i,:])**2)

        residue[j] = sum(leastsq)
        j = j + 1
    
    best_index = numpy.argmin(residue)
    best_delG = g_space[best_index]

    a = max([0, best_index-100])
    b = min([len(g_space), best_index+100])
    
    fig_fitting = plt.figure()
    plt.plot(g_space[a:b], residue[a:b])
    plt.xlabel("Free energy of sticky end (kcal/mol)")
    plt.ylabel("Sum of squared residues")
    outputobj = StringIO.StringIO()
    fig_fitting.savefig(outputobj, format="png")
    data_fitting = outputobj.getvalue().encode('base64').replace('\n', '')
    outputobj.close()
    fitting_image_tag = '<img src="data:image/png;base64,{0}" height="400" width="500">'.format(data_fitting)
    
    return best_delG, fitting_image_tag

def format_input(ic_s, om_s, t, bp1, bp2):
    ic = numpy.zeros([t, 2])
    om = numpy.zeros([t, 3])
    omr = numpy.zeros([t, 3])
    bp = numpy.zeros(2)

    for i in range(t):
        ic[i,0] = float(ic_s[i][0])
        ic[i,1] = float(ic_s[i][1])
        om[i,0] = float(om_s[i][0])
        om[i,1] = float(om_s[i][1])
        om[i,2] = float(om_s[i][2])
        omr[i,:] = om[i,:]/om[i,2]

    bp[0], bp[1] = float(bp1), float(bp2)
    return ic, omr, bp

def delG_format_input(g_custom, g_low_s, g_high_s, oh_s):
    if g_custom:
        g_low, g_high = float(g_low_s), float(g_high_s)
    elif oh_s == "2":
        g_low, g_high = -12.81, -6.88    #kcal/mol
    elif oh_s == "3":
        g_low, g_high = -14.65, -7.76
    elif oh_s == "4":
        g_low, g_high = -17.22, -8.34
    elif oh_s == "5":
        g_low, g_high = -19.06, -9.22
    return g_low, g_high

def error_output(opt):
    print "Error:"
    if opt == 0:
        print "Base pairs information missing or invalid."
    elif opt == 1:
        print "No experimental information provided."
    elif opt == 2:
        print "Invalid experimental conditions."
    elif opt == 3:
        print "Invalid free energy information."
    
def fitting_output(delG, image_tag):
    print "The free energy of sticky end (best fit) is <b>%.2f kcal/mol.</b><br/>" % delG
    print image_tag


def input_check():
    form = cgi.FieldStorage()
    g_custom = False
    inputok = False
    error_code = 0
    delG = 0
    image_tag = ""
    if ("bp1" not in form or "bp2" not in form) or not (
        is_integer(form.getfirst("bp1")) and is_integer(form.getfirst("bp2"))):
        error_code = 0
    elif "expt" not in form:
        error_code = 1
    else:
        trial = form.getlist("expt")
        trial_num = len(trial)
        ic_s = [["foo" for i in range(2)] for j in range(trial_num)]    #python list, not numpy array
        om_s = [["foo" for i in range(3)] for j in range(trial_num)]    #python list, not numpy array
        inputok = True
        for i in range(trial_num):
            j = i + 1
            ic_s[i][0] = form.getfirst("ic_dna1_%d" % j, "bad")
            ic_s[i][1] = form.getfirst("ic_dna2_%d" % j, "bad")
            om_s[i][0] = form.getfirst("om_dna1_%d" % j, "bad")
            om_s[i][1] = form.getfirst("om_dna2_%d" % j, "bad")
            om_s[i][2] = form.getfirst("om_dna12_%d" %j, "bad")
            if not (is_float(ic_s[i][0]) and is_float(ic_s[i][1]) and is_float(
                om_s[i][0]) and is_float(om_s[i][1]) and is_float(om_s[i][2])):
                error_code = 2
                inputok = False
                break
            
    if inputok and ("g_low" in form or "g_high" in form):
        if ("g_low" in form and "g_high" not in form) or ("g_low" not in form
                and "g_high" in form) or not (is_float(form.getfirst("g_low"))
                        and is_float(form.getfirst("g_high"))):
            error_code = 3
            inputok = False
        elif float(form.getfirst("g_low")) >= float(form.getfirst("g_high")):
            error_code = 3
            inputok = False
        else:
            g_custom = True
        
    if inputok:
        ic, omr, bp = format_input(ic_s, om_s, trial_num,
                                       form.getfirst("bp1"), form.getfirst("bp2"))
        g_low, g_high = delG_format_input(g_custom, form.getfirst("g_low"),
                                          form.getfirst("g_high"), form.getfirst("overhang"))
        delG, image_tag = fitting_delG(trial_num, ic, omr, bp, g_low, g_high)

    return inputok, error_code, delG, image_tag

def output_html(inputok, error_code, delG, image_tag):
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
            <li><a href="../seas/about.html">Overview</a></li>
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
        <h3>By Experiments</h3>
        
        <p class="justify">
        """
    if not inputok:
        error_output(error_code)
    else:
        fitting_output(delG, image_tag)
    print """
        </p>

        </div>

        <hr>

        <p class="name">Last modified by Daphne Ng on Jan 4th, 2014</p>

        </div>
        </body>
        </html>
        """

inputok, error_code, delG, image_tag = input_check()
output_html(inputok, error_code, delG, image_tag)                
