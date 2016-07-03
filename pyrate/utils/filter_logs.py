'''
* filter log files (should rename to filter_log.py)
* trying to use python string processing
'''

import getopt
import sys
import collections

from copy import deepcopy

'''
opts, args = getopt.getopt(sys.argv[1:], 'd:i:j:v:', ['od=','of='])
# usage checking
if opts == [] and args != []:
    print 'read usage'
    sys.exit(0)
# gettting command line options
root_direct = False
ign_list = []
jst_list = []
out_direct = False
out_fn = False
verbosity = 2       # default verbosity (don't tell about stripping rows)
# assiging and checking valid options
for opt in opts:
    if opt[0] == '-d':
        # check if is a valid directory and ensure ends in '/' or '\'
        root_direct = opt[1]
    if opt[0] == '-i':
        ign_list = opt[1]
    if opt[0] == '-j':
        jst_list = opt[1]
    if opt[0] == '--od':
        out_direct = opt[1]
    if opt[0] == '--of':
        out_fn = opt[1]
    if opt[0] == '-v':
        verbosity = int(opt[1])
# checking for valid output options
if out_direct != False and out_fn != False:
    print 'can only have one or the other or none of --od and --of'
    sys.exit(0)
# check verbosity options are correct
if not (verbosity in [1, 2, 3]):
    print 'verbositiy can only be 1 or 2. default is 1 (most verbose)'
    sys.exit(0)
'''

opts, args = getopt.getopt(sys.argv[1:], 'i:')
for opt in opts:
    if opt[0] == '-i':
        # check if is a valid directory and ensure ends in '/' or '\'
        in_fn = opt[1]

in_fp  = open(in_fn, 'r')
#out_fp = open('/home/gap/filtered.log', 'w')

'''
p_lt = 50       # filter out everything less than this percentage...

line_buffer = collections.deque(['']*10, 10)    # buffer the last 10 lines
#print line_buffer
#while True: pass

for line in in_fp:
    if line[0] == '*' or line[4] == '*':
        print line

    if line[-2] == '%':
        #print line.split(' ')[-1]
        if float(line.split(' ')[-1][:-2]) > p_lt:
            # we wan't this line and its parents...
            #print line
            pass

    # add this line to the queue

in_fp.close()
out_fp.close()
'''

#my_filter = 50  # filter on 50%

# todo: check if in file is too big to store in memory
in_lines = in_fp.readlines()
#print in_lines
in_fp.close()

# get preview in TKinter
from Tkinter import *
class FilterGUI:
    def showGuiData(self):
        # prints data in self.guiData to text area
        for line_n, line in enumerate(self.guiData):
            # first row/cols = 1/0
            pos_str = str(line_n+1)+'.0'
            self.logTxt.insert(pos_str, line)

    def clearMsgBox(self):
        self.messageBox.delete('0','end')
        self.messageBox.insert('0', '')

    def reset(self):
        # resets fields and textbox back to original
        '''
        self.tolText.delete('0','end')
        self.nanText.delete('0','end')
        '''
        self.totText.delete('0','end')
        self.guiData = deepcopy(in_lines)
        self.logTxt.delete('1.0', 'end')
        self.showGuiData()

    def apply(self):
        # apply the filters
        '''
        tol_fail = False
        nan_fail = False
        try:
            tol_filter = float(self.tolText.get())
            self.clearMsgBox()
        except ValueError:
            tol_fail = True
        try:
            nan_filter = float(self.nanText.get())
            self.clearMsgBox()
        except ValueError:
            nan_fail = True
        '''
        fail = False
        try:
            tot_filter = float(self.totText.get())
            self.clearMsgBox()
        except ValueError:
            fail = True
        if fail:
            self.messageBox.delete('0','end')
            self.messageBox.insert('0', 'invalid float encountered')
            self.totText.delete('0','end')
            return
        else:
            self.clearMsgBox()

        # got valid filters. filter current guiData
        filtered = []
        current_number = None
        for num, line in enumerate(self.guiData):
            if line[0] == '*' or line[1] == '*':
                filtered.append(line)
            if '* total fail' in line:
                # check it tolerance satisfied...
                if float(line.split(' ')[-1][:-2]) > tot_filter:

                    # iterate backwards until reach number
                    it = num-1
                    stop = False
                    while True:
                        if self.guiData[it][3] == '*' and (not stop):
                            app_from = it
                            stop = True
                        if self.guiData[it][2] == '*':
                            # this is the number line
                            this_num = self.guiData[it].split(' ')[-1]
                            break
                        it -= 1

                    if current_number != this_num:
                        current_number = this_num
                        filtered.append(self.guiData[it])
                    # append this non number stuff
                    print range(app_from, num)
                    for app_it in range(app_from, num+1):
                        filtered.append(self.guiData[app_it])

        self.guiData = filtered
        #self.guiData = ['nothing\n', 'nothing\n', 'make sammiches\n']
        self.logTxt.delete('1.0', 'end')
        self.showGuiData()

    def __init__(self):
        self.tk = Tk()
        # create frames want to pack into. todo: not sure if this is how frames and positioning work. leave till later
        # tk frames
        self.showFrame = Frame(self.tk)
        self.showFrame.pack(side='bottom')
        self.topFrame = Frame(self.tk)
        self.topFrame.pack()
        # top frames

        self.filterFrame = Frame(self.topFrame)
        self.filterFrame.pack(side='left', anchor='w')
        '''
        self.tolFilterFrame = Frame(self.filterFrame)
        self.tolFilterFrame.pack()
        self.nanFilterFrame = Frame(self.filterFrame)
        self.nanFilterFrame.pack()
        '''
        self.totFilterFrame = Frame(self.filterFrame)
        self.totFilterFrame.pack()
        # ----------------------
        self.buttonsFrame = Frame(self.topFrame)
        self.buttonsFrame.pack()
        # ----------------------
        self.messageFrame = Frame(self.buttonsFrame)
        self.messageFrame.pack()

        # add tolerance filter label and text box
        '''
        self.tolLabel = Label(self.tolFilterFrame, text='tol filter:')
        self.tolLabel.pack(side='left', anchor='w')
        # -----------------------
        self.tolText = Entry(self.tolFilterFrame, width='25')
        self.tolText.pack(anchor='w')
        # =======================
        self.nanLabel = Label(self.nanFilterFrame, text='NaN filter:')
        self.nanLabel.pack(side='left', anchor='w')
        # -----------------------
        self.nanText = Entry(self.nanFilterFrame, width='25')
        self.nanText.pack(anchor='w')
        # =======================
        '''
        self.totLabel = Label(self.totFilterFrame, text='tot filter:')
        self.totLabel.pack(side='left', anchor='w')
        # -----------------------
        self.totText = Entry(self.totFilterFrame, width='25')
        self.totText.pack(anchor='w')


        # add buttons
        self.resetButton = Button(self.buttonsFrame, text='reset', command=self.reset)
        self.resetButton.pack(anchor='w')
        self.applyButton = Button(self.buttonsFrame, text='apply', command=self.apply)
        self.applyButton.pack(anchor='w')
        self.messageLabel = Label(self.buttonsFrame, text='Message:')
        self.messageLabel.pack(side='left')
        # -----------------------
        self.messageBox = Entry(self.buttonsFrame, width='30')
        self.messageBox.pack()

        # show the config file in text area...
        self.logTxt = Text(self.showFrame, width=100, height=50)
        self.logTxt.pack()

        # actually show it
        self.guiData = deepcopy(in_lines)
        self.showGuiData()

        # kick it off
        self.tk.mainloop()

gui = FilterGUI()