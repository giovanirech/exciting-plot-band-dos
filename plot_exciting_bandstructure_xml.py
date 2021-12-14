#!/usr/bin/env python3
# plot_exciting_bandstructure_xml.py
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse

def read_bands_Ldiff(xml_file):
    tree = ET.parse(xml_file)
    bandstrucutre = tree.getroot()
    
    #Vamos criar um dicionário com as bandas separadas por elemento (species)
    
    dictionary_species = {}
    for species_i in range(len(bandstrucutre.findall('species'))):
        
        key_element_symbol = bandstrucutre.findall('species')[species_i].attrib['chemicalSymbol']
    
        bands = bandstrucutre.findall('species')[species_i].getchildren()[0].getchildren() #todas as bandas do atomo
        #       ^root         ^species                      ^atom            ^bands
        
        
        
        #lista de dicionarios. Cada dicionario é uma banda
        lista_bandas = []
        for band in bands:
            distance = np.array([float(i.get('distance')) for i in band.getchildren()])
            energy = np.array([float(i.get('eval')) for i in band.getchildren()])
            suml = np.array([float(i.get('sum')) for i in band.getchildren()])
            inters = np.array([1-float(i.get('sum')) for i in band.getchildren()])
            s = np.array([float(i.getchildren()[0].get('character')) for i in band.getchildren()])
            p = np.array([float(i.getchildren()[1].get('character')) for i in band.getchildren()])
            d = np.array([float(i.getchildren()[2].get('character')) for i in band.getchildren()])
            f = np.array([float(i.getchildren()[3].get('character')) for i in band.getchildren()])
            g = np.array([float(i.getchildren()[4].get('character')) for i in band.getchildren()])
            banda = {'distance': np.array(distance),
                     'energy'  : np.array(energy),
                     'sum'     : np.array(suml),
                     'inters'  : np.array(inters),
                     's'       : np.array(s),
                     'p'       : np.array(p),
                     'd'       : np.array(d),
                     'f'       : np.array(f),
                     'g'       : np.array(g)}
            lista_bandas.append(banda)
        
        dictionary_species[key_element_symbol] = lista_bandas
    ##############################################
    # Outros elementos gráficos
    vertice_labels = []
    vertice_coords = []
    for vertice in bandstrucutre.findall('vertex'):
        vertice_coords.append(float(vertice.get('distance')))
        label = vertice.get('label')
        if label == 'GAMMA':
            label = '$\Gamma$'
        vertice_labels.append(label)
    return dictionary_species, vertice_coords, vertice_labels

#from matplotlib.colors import to_rgba_array, to_rgb

def plot_bandstructure_atoms(dictionary_species, vertice_coords, vertice_labels, 
                             l_resolved ,output_file,elims, efermi, orbthickness, orbalpha, orbcolors,
                             linethickness, linealpha, linecolor):
    
    #plt.rc('font', family='QuattroCento Sans')
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)
    plt.rc('axes', labelsize=20)
    #plt.rc('mathtext', rm='Quattrocento Sans',it='Quattrocento Sans:italic', bf='Quattrocento Sans:bold', fontset='custom' )
    
    elements = dictionary_species.keys()
    N = len(elements)
    
    
    fig = plt.figure(figsize=(10,4*N), dpi=150)
    
    for n,element in enumerate(elements):
        #PLOT FIRST VOLUME
        ax0 = fig.add_subplot(N,1,n+1)
    
        lista_bandas = dictionary_species[element]
    
        #get x limits from data
        xmin,xmax = lista_bandas[0]['distance'][0],lista_bandas[0]['distance'][-1]
    
        ##### Diferenciando L
        #draw curves
        if l_resolved:
            for banda in lista_bandas:
                alpha_geral = orbalpha
                multiplier = orbthickness #define a espessura da banda de acordo com o orbital
                ax0.plot(banda['distance'], banda['energy'], color='grey', lw=linethickness, alpha=linealpha)
                orbitals = ['s',  'p',  'd',  'f',  'g']
                colors =   orbcolors
                for i in range(len(orbitals)):
                    l = ax0.fill_between(banda['distance'],
                                         banda['energy']+multiplier*banda[orbitals[i]],
                                         banda['energy']-multiplier*banda[orbitals[i]],
                                         color=colors[i], alpha = alpha_geral, lw=0.0)
        
            #dummy plot for legend purposes
            for i in range(len(orbitals)):
                    l = ax0.fill_between([-1], [1],[0], color=colors[i], alpha = alpha_geral, label=orbitals[i])
            ax0.legend(loc=1, ncol=5, fontsize=12)
        
        else: #SEM DIFERENCA ENTRE L
            for banda in lista_bandas:
                ax0.plot(banda['distance'], banda['energy'], color=linecolor, lw=linethickness, alpha=linealpha)
                
                
        ################################
        #get y limits for later
        ymin, ymax = elims
        if ymin==None:
            ymin = plt.ylim()[0]
        if ymax==None:
            ymax = plt.ylim()[1]
        
        #draw fermi level
        ax0.plot([xmin,xmax],[efermi,efermi], color='k', alpha = 0.5, linestyle = '--')
        #ax1.annotate('E$_\mathrm{F}$', xy = (0,ymax-4), fontsize=20)
        ax0.annotate(str(element), 
                 xycoords = 'axes fraction', xy=(0.02, 0.95),
                 fontsize=16,
                 horizontalalignment='left',
                 verticalalignment='top')
        #draw high symmetry lines
        ax0.vlines(vertice_coords, ymin, ymax, color='k', linewidth = 0.2)
        #adjust limits
        ax0.set_xlim(xmin, xmax)
        ax0.set_ylim(ymin, ymax)
        ax0.set_xticks([])
    
    
    #set positions and labels of x axis
    plt.xticks(vertice_coords, labels = vertice_labels)
    fig.text(0.02, 0.5, 'Energy (Ha)', va='center', rotation='vertical', fontsize = 20)
    plt.subplots_adjust(hspace=.02)
    #plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')
    
#========================================================================================    
if __name__=='__main__':
    # Create the parser
    parser = argparse.ArgumentParser(description = 'Plot the electronic band structures from a XML file generated by the Exciting code. The program creates one panel for each atom found in the input file showing the electronic band structure projected upon that atom. By default, the program resolves the contribution for each type of orbital (s,p,d,f,g) and how much of the electronic density is of that type at each k-vector. This in encoded in the band structure plot trhough the lines thickness and colors. Each color represent a orbital type while its thickness express the relative character of that type. For example, a band at a specific k-vector with strong d (green) character and weak s (blue) character will be expressed in the plot as a line composed of a thick green band and a thin blue band overlayered. Notice that the interstitial contribution is NOT represented in the bandstructure plot and, therefore, a thin line indicates that, at that k-point, the electronic contribution is mostly interstitial. The thickness and opacity of such bands can be controlled using the -t and -a options, respectively.  ', epilog = '--Created by G.L.Rech')
    # Add the arguments
    parser.add_argument('xml_file', metavar='xml_file', type=str,
                       help='Path to the xml file.')
    parser.add_argument('-N','--Noorbital', action = 'store_false',
                       help = 'Supress the diferentiation by orbitals. The default is to differentiate by s,p,d,f, and g orbitals')
    parser.add_argument('-o','--output', type=str,
                       help='Name of the output file. The extension can be used to change the type of the output file. The default is <input_file>.pdf')
    parser.add_argument('--Emin', type=float, action='store',
                       help='The minimum value of energy in the plot (in Ha). The default is automatically chosen to show all the data available.')
    parser.add_argument('--Emax', type=float, action='store',
                       help='The maximum value of energy in the plot (in Ha). The default is automatically chosen to show all the data available.')
    parser.add_argument('-F', '--fermi', type=float, action='store', default=0.0,
                       help='The fermi energy, in Ha. The horizontal dashed line will be placed at this value. Default is 0.0.')
    parser.add_argument('-T','--orbthickness', type=float, action='store', default=0.1,
                       help='This value controls the relative thickness of the bands when L-resolved. the default is 0.1.')
    parser.add_argument('-t','--linethickness', type=float, action='store', default=0.2,
                       help='Controls the thickness of the main lines of the bandstructure. When the -N option is used, this is the only line shown. The default value is 0.2')
    parser.add_argument('-A','--orbalpha', type=float, action='store', default=0.3,
                       help='This value controls the opacity of the bands that characterizes the color coded orbitals. The default is 0.3. It should be in the range [0.0,1.0].')
    parser.add_argument('-a','--linealpha', type=float, action='store', default=0.5,
                       help='This value controls the opacity of the main lines of the bandstructure. The default is 0.5. It should be in the range [0.0,1.0].')
    parser.add_argument('-C','--orbcolors', nargs=5, default = ['C0', 'C1', 'C2', 'lightgray', 'C4'],
                       help='Define the colors of the s, p, d, f, and g bands. Notice that, if used, 5 colors should be passed as arguments. This list can be any color format accepted by matplotlib. The default colors are [C0, C1, C2, lightgray, C4].')
    
    parser.add_argument('-c','--linecolor', default = 'k',
                       help='The color of the bandstructure lines. This is the color of the main lines in the plot. If the -N option is used, this is the only color shown. Default is black')
    
    
    #get the passed arguments
    args = parser.parse_args()
    
    if args.output == None:
        args.output=args.xml_file.rstrip('.'+args.xml_file.split('.')[-1])+'.pdf'

    dictionary_species, vertice_coords, vertice_labels = read_bands_Ldiff(args.xml_file)
    plot_bandstructure_atoms(dictionary_species, vertice_coords, vertice_labels,
                             l_resolved = args.Noorbital,
                             output_file = args.output,
                             elims = (args.Emax, args.Emin),
                             efermi = args.fermi,
                             orbthickness = args.orbthickness,
                             orbalpha = args.orbalpha, 
                             orbcolors = args.orbcolors,
                             linethickness = args.linethickness,
                             linealpha = args.linealpha,
                             linecolor = args.linecolor)
   