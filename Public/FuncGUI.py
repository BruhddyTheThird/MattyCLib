import tkinter as tk
from tkinter import ttk, messagebox
import sympy as sp
import os, sys
# go up one directory level from this file's directory:
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# prepend parent directory to the system path:
sys.path.insert(0, path)
import FunctionLib.MTHFunctions as mtf

root = tk.Tk(None,baseName=None,className='JoemApp',useTk=1)
Label1 = tk.Label(root,text='Function List:')
Label1.pack()

tk.mainloop()