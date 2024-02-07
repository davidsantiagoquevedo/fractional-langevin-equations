# -*- coding: utf-8 -*-
"""
Created on Tue Feb 7 2024
@author: davidsantiagoquevedo
"""

def add_caption_letter(ax, cap):
    ax.text(.02, .98, cap, ha='left', va='top', transform=ax.transAxes)