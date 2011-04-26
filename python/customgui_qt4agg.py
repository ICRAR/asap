import os
import matplotlib, numpy
from asap.logging import asaplog, asaplog_post_dec
from matplotlib.patches import Rectangle
from asap.parameters import rcParams
from asap._asap import stmath
from asap.customgui_base import *

import PyQt4 as qt

class CustomToolbarQT4Agg(CustomToolbarCommon,  qt.QtGui.QToolBar):
    def __init__(self,parent):
        from asap.asapplotter import asapplotter
        if not isinstance(parent,asapplotter):
            return False
        if not parent._plotter:
            return False
        self._p = parent._plotter
        self.figmgr = self._p.figmgr
        self.canvas = self.figmgr.canvas
        self.mode = ''
        self.button = True
        self.pagecount = None
        CustomToolbarCommon.__init__(self,parent)
#         self.notewin = NotationWindowQT4Agg(master=self.canvas)
        self._add_custom_toolbar()

    def _add_custom_toolbar(self):
        qt.QtGui.QToolBar.__init__(self,parent=self.figmgr.window)
        self.figmgr.window.addToolBar(qt.QtCore.Qt.BottomToolBarArea,self)
        self.bNote = self._NewButton(master=self,
                                     text='notation',
                                     command=self.modify_note)
        self.bNote.setCheckable(True)

        self.bStat = self._NewButton(master=self,
                                     text='statistics',
                                     command=self.stat_cal)
        self.bStat.setCheckable(True)

        # page change oparations
        frPage = qt.QtGui.QWidget(parent=self,flags=qt.QtCore.Qt.Tool)
        loPage = qt.QtGui.QHBoxLayout(self)
        loPage.addStretch(1)
        self.lPagetitle = qt.QtGui.QLabel('Page:',parent=frPage)
        self.lPagetitle.setMargin(5)
        loPage.addWidget(self.lPagetitle)
        self.pagecount = qt.QtGui.QLabel(parent=frPage)
        self.pagecount.setStyleSheet("background-color: white")
        self.pagecount.setMargin(3)
        self.pagecount.setText('   1')
        loPage.addWidget(self.pagecount)
        
        self.bNext = self._NewButton(master=frPage,
                                     text=' + ',
                                     command=self.next_page,addparent=False)
        loPage.addWidget(self.bNext)
        self.bPrev = self._NewButton(master=frPage,
                                     text=' - ',
                                     command=self.prev_page,addparent=False)
        loPage.addWidget(self.bPrev)
        frPage.setLayout(loPage)
        self.addWidget(frPage)

        self.bQuit = self._NewButton(master=self,
                                     text='Quit',
                                     command=self.quit)

#         if os.uname()[0] != 'Darwin':
#             self.bPrev.config(padx=5)
#             self.bNext.config(padx=5)

        self.pagecount.setText(' '*4)

        self.disable_button()
        return

    def _NewButton(self, master, text, command,addparent=True):
        b = qt.QtGui.QPushButton(text,parent=master)
        if addparent: master.addWidget(b)
        master.connect(b,qt.QtCore.SIGNAL('clicked()'),command)
#         if os.uname()[0] == 'Darwin':
#             b = Tk.Button(master=master, text=text, command=command)
#         else:
#             b = Tk.Button(master=master, text=text, padx=2, pady=2,
#                           command=command)
        return b

    def show_pagenum(self,pagenum,formatstr):
        self.pagecount.setText(formatstr % (pagenum))

    def spec_show(self):
        if not self.figmgr.toolbar.mode == '' or not self.button: return
        self.figmgr.toolbar.set_message('spec value: drag on a spec')
        if self.mode == 'spec': return
        self.mode = 'spec'
#         self.notewin.close_widgets()
        self.__disconnect_event()
        self._p.register('button_press',self._select_spectrum)

    def stat_cal(self):
        if not self.figmgr.toolbar.mode == '' or not self.button:
            # Get back button status BEFORE clicked
            self.bStat.setChecked(not self.bStat.isChecked())
            return
        if self.mode == 'stat':
            # go back to spec mode
            self.bStat.setChecked(False)
            self.spec_show()
            return
        self.figmgr.toolbar.set_message('statistics: select a region')
        self.bStat.setChecked(True)
        self.bNote.setChecked(False)
        self.mode = 'stat'
#         self.notewin.close_widgets()
        self.__disconnect_event()
        self._p.register('button_press',self._single_mask)

    def modify_note(self):
        if not self.figmgr.toolbar.mode == '':
            # Get back button status BEFORE clicked
            self.bNote.setChecked(not self.bNote.isChecked())
            return
        self.figmgr.toolbar.set_message('text: select a position/text')
        if self.mode == 'note':
            self.bNote.setChecked(False)
            self.mode = 'none'
            self.spec_show()
            return
        self.bStat.setChecked(False)
        self.bNote.setChecked(True)
        self.mode = 'note'
        self.__disconnect_event()
        self._p.register('button_press',self._mod_note)

    def quit(self):
        self.__disconnect_event()
        self.disable_button()
        self._p.quit()

    def enable_button(self):
        if self.button: return
        self.bStat.setEnabled(True)
        self.button = True
        self.spec_show()

    def disable_button(self):
        if not self.button: return
        self.bStat.setChecked(False)
        self.bStat.setDisabled(True)
        self.button = False
        self.mode = ''
        self.__disconnect_event()

    def enable_next(self):
        self.bNext.setEnabled(True)

    def disable_next(self):
        self.bNext.setDisabled(True)

    def enable_prev(self):
        self.bPrev.setEnabled(True)

    def disable_prev(self):
        self.bPrev.setDisabled(True)

    # pause buttons for slow operations
    def _pause_buttons(self,operation="end",msg=""):
        buttons = ["bStat","bNote","bQuit"]
        if operation == "start":
            enable = False
        else:
            enable = True
        for btn in buttons:
            getattr(self,btn).setEnabled(enable)
        self.figmgr.toolbar.set_message(msg)

    def delete_bar(self):
        self.__disconnect_event()
        self.destroy()

    def __disconnect_event(self):
        self._p.register('button_press',None)
        self._p.register('button_release',None)


