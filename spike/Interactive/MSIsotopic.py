#!/usr/bin/env python 
# encoding: utf-8

"""
A utility to use the isotope module with spike within jupyter


First version MAD July 2021
"""

from ipywidgets.widgets.widget_box import VBox
import matplotlib.pylab as plt
import numpy as np

from ipywidgets import Text, Layout, Output, Box, Label, Button, HTML, BoundedIntText, Dropdown, FloatText
import ipywidgets as widgets
from IPython.display import display, Javascript, Markdown
from traitlets.traitlets import Float

from spike.Interactive import INTER as I
try:
    from isotope import isotopes as iso
    isavailable = True
except ImportError:
    isavailable = False

if not isavailable:
    print("""
*** MSIsotopic not available

The MSIsotopic module provides a set of tools for computing MS spectra from atomic formula and peptide sequences.
You can do without it, if you want to have to get it
- at https://github.com/delsuc/isotope from the standard operations
- an additional tool for fine isotopic structures uses *neutronstar* from P. Kreitzberg et al (Montana Univ.)
  available at: https://bitbucket.org/orserang/neutronstar.git

""")
    raise ImportError


iso.NS = "./isotope/neutronstar/neutronstar"


class SpIsotope(I.Show1D):
    def __init__(self, data):
        super().__init__( data, title=None, figsize=None, reverse_scroll=False, show=True, create_children=True)

def fspacer(n):
    "used in flex"
    return Label(' ',layout=Layout(flex='%d'%(n), width='auto'))
def ftext(text,n):
    "used in flex"
    return Label(text,tooltip='coucou',layout=Layout(flex='%d'%(n), width='auto'))

class Isotope():
    """
    this class does the whole thing !

    provides 3 widgets:
    entrybar: to enter the molecule definition
    resbar: to list the isotopic pattern
    drawbar: to draw the simulated spectra

    and a completebar to rule them all

    finally .full() implements the complete display
    """
    def __init__(self, axis=None):
        "axis is an optional matplotlib axis on which to draw the spectra"

        self.axis = axis
        # Entry bar
        self.meth = Dropdown( options=['atomic formula', '1 letter peptide'],
            value='atomic formula',
            layout=Layout(flex='2', width='auto'))
        self.formula = Text( value='',
            placeholder='(CH)5 C CH3 ',
            continuous_update=False,
            layout=Layout(flex='15', width='auto'))
        self.formula.observe(self.compute)
        self.meth.observe(self.default_form)
        self.quality = Dropdown( options=['isotopic', 'fine isotopic'],
            value='isotopic',
            layout=Layout(flex='2', width='auto'))
        self.quality.observe(self.compute)
        self.docinput = Button(description='?',
            tooltip="documentation on entering formula",
            layout=Layout(flex='1', width='auto'))
        self.docinput.on_click(self.showdoc)
        # flex 20

        # draw bar
        self.todraw = Dropdown( options=['bar', 'spectrum', 'both'], value='bar',
            layout=Layout(flex='3', width='auto'))
        self.todraw.observe(self.selectRP)
        self.charge = BoundedIntText( value=1,
            description='charge',
            min=-100, max=100, step=1,
            layout=Layout(flex='1', width='10em'))
        self.RP = BoundedIntText( value=30,
            min=1, max=10000, step=10,
            description='RP',
            tooltip='Resolving Power (in thousand)',
            disabled = True,
            layout=Layout(flex='1', width='12em'))
        self.scale = FloatText(value=1E5,
            description='scaling',
            layout=Layout(flex='1', width='14em'))
        self.color = Dropdown( options=I.Colors, value='red',
            layout=Layout(flex='2', width='auto'))
        self.bdraw = Button(description='Draw', button_style='success',
            tooltip="Draw simulated isotopic patterns",
            layout=Layout(flex='2', width='auto'))
        self.bdraw.on_click(self.draw)
        # flex 10

        # result bar
        self.atom_form = HTML(
            value='<span style="color: grey;font-style: italic;">results will come here</span>',
            layout=Layout(flex='10', width='auto'))
        self.bmass = Button(description='Sorted by Mass',
            tooltip="List isotopes by increasing mass value",
            layout=Layout(flex='2', width='auto'))
        self.bmass.on_click(self.masstable)
        self.babund = Button(description='Sorted by Abundance',
            tooltip="List isotopes by decreasing abundance",
            layout=Layout(flex='2', width='auto'))
        self.babund.on_click(self.abundtable)
        self.table = Output(width='100%')
        self.bcleartable = Button(description='clear table',
            tooltip="Clear the on-screen list",
            disabled=True,
            layout=Layout(flex='2', width='auto'))
        self.bcleartable.on_click(self.cleartable)
        #flex 16

        hbox_layout = Layout(display='flex',
                    flex_flow='row',
                    align_items='stretch',
                    width='90%')
        vbox_layout = Layout(display='flex',
                    flex_flow='column',
                    align_items='stretch',
                    width='90%')
        self.entrybar = Box([ftext("Enter formula:",3), self.meth, self.formula, self.docinput, self.quality],
                         layout=hbox_layout)
        self.drawbar = Box([fspacer(5),  ftext("Draw isotopic pattern as",8), self.todraw,
                self.scale, self.charge, self.RP, ftext("k",1), 
                self.color, self.bdraw],
                # ftext("charge",2), self.charge,
                # ftext("RP",2), self.RP, ftext("k",1), 
                # ftext("scaling",2), self.scale, 
                # self.color, self.bdraw],
                        layout=hbox_layout)
        self.resbar = Box(  [Box([self.atom_form, self.bmass, self.babund, self.bcleartable], layout=hbox_layout),
                             self.table],
                            layout=vbox_layout)
        self.completebar = Box([self.entrybar, self.drawbar, self.resbar], layout=vbox_layout)

    def full(self, data):
        "a composite view"
        if self.axis is not None:
            raise Exception('You should not provide an axis for this feature to work')
        Show = I.Show1D(data)
        self.axis = Show.ax
        vbox_layout = Layout(display='flex',
            flex_flow='column',
            align_items='stretch',
            width='100%')
        self.completebar = Box([self.entrybar, self.drawbar, Show.fig.canvas, self.resbar], layout=vbox_layout)
        return self.completebar
    def showdoc(self, e):
        if self.meth.value == 'atomic formula':
            text = iso.parse_formula.__doc__
        else:
            text = iso.parse_peptide.__doc__
        I.jsalert("\\n".join( text.split('\n') ))
    def default_form(self, e):
        "change placeholder text for formula entry"
        if self.meth.value == 'atomic formula':
            self.formula.placeholder = 'eg: CH3 CH2 OH'
        else:
            self.formula.placeholder = 'one letter coded PEPTIDE'
    def selectRP(self, e):
        if self.todraw.value == 'bar':
            self.RP.disabled = True
        else:
            self.RP.disabled = False
    def draw(self, e=None):
        ""
        self.compute()
        if self.axis is None:
            self.fig, self.axis = plt.subplots()
        if self.todraw.value in ('bar', 'both'):
            self.Distrib.bar(largest=self.scale.value, charge=abs(self.charge.value),
                color=self.color.value,)
        if self.todraw.value in ('spectrum', 'both'):
            self.Distrib.draw(RP=1000*self.RP.value, largest=self.scale.value, charge=abs(self.charge.value), 
                color=self.color.value, label="simulated pattern - RP=%dk"%self.RP.value)

    def compute(self, e=None):
        "computes every thing - called each time something is changed in entrybar"
        if e is not None:
            if e['name'] != 'value':
                return
        try:
            if self.meth.value == 'atomic formula':
                self.molecule = iso.parse_formula(self.formula.value)
            else:
                self.molecule = iso.parse_peptide(self.formula.value, extended=True)
        except:
            self.meth.value
            I.jsalert("%s not a valid %s\\nError in Molecular Definition\\nCheck documentation ( [?] button) "% \
                    (self.formula.value, self.meth.value))
            return
        self.atom_form.value = "%s  monoisotopic mass: %.6f   average mass: %.3f"% \
                                       (iso.HTMLformula(self.molecule), self.molecule.monoisotop(), self.molecule.average())
        if self.molecule.average() >0:  # would be empty string
            if self.quality.value == 'fine isotopic':
                try:
                    self.Distrib = self.molecule.fineisotopicdistrib()
                except Exception as e:
                    I.jsalert(str(e))
            else:
                self.Distrib = self.molecule.distribution()
        self.table.clear_output()

    def masstable(self, e):
        "list by masses"
        self.Distrib.sort_by_mass()
        self.table.clear_output()
        with self.table:
            print(self.Distrib)
        self.bcleartable.disabled = False
    def abundtable(self, e):
        "list by abundance"
        self.Distrib.sort_by_intens()
        self.table.clear_output()
        with self.table:
            print(self.Distrib)
        self.bcleartable.disabled = False
    def cleartable(self, e):
        "clear table"
        self.table.clear_output()
        self.bcleartable.disabled = True
