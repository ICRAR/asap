import os
import matplotlib
from asap.logging import asaplog, asaplog_post_dec

######################################
##    Notation box window           ##
######################################
class NotationWindowCommon:
    def __init__(self,master=None):
        #self.parent = master
        self.canvas = master
        self.event = None
        self.note = None
        self.ancval = None
        self.anchors = ["figure","axes","data"]
        self.seltext = {}
        self.numnote = 0

    @asaplog_post_dec
    def print_text(self):
        anchor = self.anchors[self.ancval.get()]
        notestr = self._get_note().rstrip("\n")
        if len(notestr.strip()) == 0:
            self._clear_textbox()
            #print "Empty string!"
            return

        myaxes = None
        calcpos = True
        xpos=None
        ypos=None
        if self.seltext:
            # You are modifying a text
            mycanvas = self.canvas
            oldanch = self.seltext['anchor']
            if oldanch != 'figure':
                myaxes = self.seltext['parent']
            calcpos = (anchor != oldanch)
            if not calcpos:
                # printing text in the same coord.
                # you don't have to recalc position
                (xpos, ypos) = self.seltext['textobj'].get_position()
                transform = self.seltext['textobj'].get_transform()
                parent = self.seltext['parent']
            elif anchor == "figure":
                # converting from "axes"/"data" -> "figure"
                (x, y) = self.seltext['textobj']._get_xy_display()
            elif oldanch == "data":
                # converting from "data" -> "axes".
                # need xdata & ydata in the axes
                (x, y) = self.seltext['textobj'].get_position()
            else:
                # converting "figure"/"axes" -> "data"
                # need to calculate xdata & ydata in the axes
                pixpos = self.seltext['textobj']._get_xy_display()
                (w,h) = mycanvas.get_width_height()
                relpos = (pixpos[0]/float(w), pixpos[1]/float(h))
                if not myaxes:
                    myaxes = self._get_axes_from_pos(relpos,mycanvas)
                    if not myaxes:
                        raise RuntimeError, "Axes resolution failed!"
                (x, y) = self._convert_pix2dat(relpos,myaxes)
            self._remove_seltext()
        elif self.event:
            mycanvas = self.event.canvas
            myaxes = self.event.inaxes
            if myaxes and (anchor != "figure"):
                x = self.event.xdata
                y = self.event.ydata
            else:
                x = self.event.x
                y = self.event.y
        else:
            raise RuntimeError, "No valid position to print data"
            return

        # now you know 
        picker = True
        # alignment of the text: ha (horizontal), va (vertical)
        ha = 'left'
        va = 'center'
        if not calcpos:
            # you aready know parent, tansform, xpos and ypos
            pass
        elif anchor == "figure":
            # text instance will be appended to mycanvas.figure.texts
            parent = mycanvas.figure
            transform = parent.transFigure
            (w,h) = mycanvas.get_width_height()
            xpos = x/float(w)
            ypos = y/float(h)            
        elif myaxes:
            ## text instance will be appended to myaxes.texts
            parent = myaxes
            if anchor == "axes":
                transform = myaxes.transAxes
                lims = myaxes.get_xlim()
                xpos = (x-lims[0])/(lims[1]-lims[0])
                lims = myaxes.get_ylim()
                ypos = (y-lims[0])/(lims[1]-lims[0])
            else:
                # anchored on "data"
                transform = myaxes.transData
                xpos = x
                ypos = y
        parent.text(xpos,ypos,notestr,transform=transform,
                    ha=ha,va=va,picker=picker)
        mycanvas.draw()

        self.numnote += 1

        self._clear_textbox()
        msg = "A note added to figure: str = '"+notestr+"'"
        msg += "   at ["+str(xpos)+", "+str(ypos)+"] in "+anchor+"-coordinate"
        msg += "\nTotal number of notes are "+str(self.numnote)
        asaplog.push( msg )

    def _get_axes_from_pos(self,pos,canvas):
        if len(pos) != 2:
            raise ValueError, "pixel position should have 2 elements"
        for axes in canvas.figure.axes:
            ##check if pos is in the axes
            #if axes.contains_point(pos): ### seems not working
            #    return axes
            axbox = axes.get_position().get_points()
            if (axbox[0][0] <= pos[0] <= axbox[1][0]) and \
               (axbox[0][1] <= pos[1] <= axbox[1][1]):
                return axes
        return None
        
    def _convert_pix2dat(self,pos,axes):
        # convert a relative position from lower-left of the canvas
        # to a data in axes
        if len(pos) != 2:
            raise ValueError, "pixel position should have 2 elements"
        # left-/bottom-pixel, and pixel width & height of the axes
        bbox = axes.get_position()
        lbpos = bbox.get_points()[0]
        wax = bbox.width
        hax = bbox.height
        # check pos value
        if (pos[0] < lbpos[0]) or (pos[1] < lbpos[1]) \
               or (pos[0] > (lbpos[0]+wax)) or (pos[1] > (lbpos[1]+hax)):
            raise ValueError, "The position is out of the axes"
        xlims = axes.get_xlim()
        ylims = axes.get_ylim()
        wdat = xlims[1] - xlims[0]
        hdat = ylims[1] - ylims[0]
        xdat = xlims[0] + wdat*(pos[0] - lbpos[0])/wax
        ydat = ylims[0] + hdat*(pos[1] - lbpos[1])/hax
        return (xdat, ydat)

    @asaplog_post_dec
    def _selected_text(self,event):
        (w,h) = event.canvas.get_width_height()
        dist2 = w*w+h*h
        selected = {}
        for textobj in self.canvas.figure.texts:
            if textobj.contains(event)[0]:
                d2 = self._get_text_dist2(event,textobj)
                if dist2 >= d2:
                    dist2=d2
                    selected={'anchor': 'figure', \
                             'parent': event.canvas.figure, 'textobj': textobj}
                    msg = "Fig loop: a text, '"+textobj.get_text()+"', at "
                    msg += str(textobj.get_position())+" detected"
                    asaplog.push(msg)
        for ax in self.canvas.figure.axes:
            for textobj in ax.texts:
                if textobj.contains(event)[0]:
                    d2 = self._get_text_dist2(event,textobj)
                    if dist2 >= d2:
                        anchor='axes'
                        if ax.transData == textobj.get_transform():
                            anchor='data'                    
                        selected={'anchor': anchor, 'parent': ax, 'textobj': textobj}
                        msg = "Ax loop: a text, '"+textobj.get_text()+"', at "
                        msg += str(textobj.get_position())+" detected"
                        asaplog.push(msg)

        return selected

    def _get_text_dist2(self,event,textobj):
        (x,y) = textobj._get_xy_display()
        return (x-event.x)**2+(y-event.y)**2

    def delete_note(self):
        #print "You selected 'OK'"
        self._remove_seltext()
        self.canvas.draw()

    @asaplog_post_dec
    def _remove_seltext(self):
        if len(self.seltext) < 3:
            raise ValueError, "Don't under stand selected text obj."
            return
        try:
            self.seltext['textobj'].remove()
        except NotImplementedError:
                self.seltext['parent'].texts.pop(self.seltext['parent'].texts.index(self.seltext['textobj']))
        self.numnote -= 1

        textobj = self.seltext['textobj']
        msg = "A note deleted from figure: str = '"+textobj.get_text()+"'"
        msg += "   at "+str(textobj.get_position())\
               +" in "+self.seltext['anchor']+"-coordinate"
        msg += "\nTotal number of notes are "+str(self.numnote)
        asaplog.push( msg )

        self.seltext = {}

    def cancel_delete(self):
        #print "You selected 'CANCEL'"
        self.seltext = {}


#####################################
##    Backend dependent Classes    ##
#####################################
### TkAgg
if matplotlib.get_backend() == 'TkAgg':
    import Tkinter as Tk
    import tkMessageBox

class NotationWindowTkAgg(NotationWindowCommon):
    def __init__(self,master=None):
        self.parent = master._tkcanvas
        NotationWindowCommon.__init__(self,master=master)
        self.textwin=self._create_textwindow(master=None)
        self.menu=self._create_modmenu(master=self.parent)
        
    def _create_textwindow(self,master=None):
        twin = Tk.Toplevel(padx=3,pady=3)
        twin.title("Annotation")
        twin.resizable(width=True,height=True)
        self.textbox = self._NotationBox(parent=twin)
        self.radio = self._AnchorRadio(parent=twin)
        self.actionbs = self._ActionButtons(parent=twin)
        
        self.textbox.pack(side=Tk.TOP,fill=Tk.BOTH,expand=True)
        self.actionbs.pack(side=Tk.BOTTOM)
        self.radio.pack(side=Tk.BOTTOM)
        #twin.deiconify()
        #twin.minsize(width=twin.winfo_width(),height=twin.winfo_height())
        twin.withdraw()
        return twin

    def _NotationBox(self,parent=None):
        textbox = Tk.Text(master=parent,background='white',
                          height=2,width=20,cursor="xterm",
                          padx=2,pady=2,undo=True,maxundo=10)
        return textbox

    def _AnchorRadio(self,parent=None):
        radio = Tk.LabelFrame(master=parent,text="anchor",
                            labelanchor="nw",padx=5,pady=3)
        self.ancval = Tk.IntVar(master=radio,value=0)
        self.rFig = self._NewRadioButton(radio,"figure",state=Tk.NORMAL,
                                         variable=self.ancval,value=0,
                                         side=Tk.LEFT)
        self.rAxis = self._NewRadioButton(radio,"panel",state=Tk.DISABLED,
                                          variable=self.ancval,value=1,
                                          side=Tk.LEFT)
        self.rData = self._NewRadioButton(radio,"data",state=Tk.DISABLED,
                                          variable=self.ancval,value=2,
                                          side=Tk.LEFT)
        # set initial selection "figure"
        self.ancval.set(0)
        return radio

    def _NewRadioButton(self,parent,text,state=Tk.NORMAL,variable=None,value=None,side=Tk.LEFT):
        rb = Tk.Radiobutton(master=parent,text=text,state=state,
                          variable=variable,value=value)
        rb.pack(side=side)
        return rb

    def _enable_radio(self):
        self.rAxis.config(state=Tk.NORMAL)
        self.rData.config(state=Tk.NORMAL)
        #self.rFig.config(state=Tk.NORMAL)
        self.rFig.select()

    def _reset_radio(self):
        self.rAxis.config(state=Tk.DISABLED)
        self.rData.config(state=Tk.DISABLED)
        self.rFig.config(state=Tk.NORMAL)
        self.rFig.select()

    def _select_radio(self,selection):
        if not selection in self.anchors:
            return
        if selection == "data":
            self.rData.select()
        elif selection == "axes":
            self.rAxis.select()
        else:
            self.rFig.select()

    def _get_note(self):
        return self.textbox.get("1.0",Tk.END)

    def _clear_textbox(self):
        self.textbox.delete("1.0",Tk.END)

    def _set_note(self,note=None):
        self._clear_textbox()
        if len(note) >0:
            self.textbox.insert("1.0",note)

    def _ActionButtons(self,parent=None):
        actbuts = Tk.Frame(master=parent)
        bCancel = self._NewButton(actbuts,"cancel",self._cancel_text,side=Tk.LEFT)
        bPrint = self._NewButton(actbuts,"print", self._print_text,side=Tk.LEFT)
        return actbuts

    def _NewButton(self, parent, text, command, side=Tk.LEFT):
        if(os.uname()[0] == 'Darwin'):
            b = Tk.Button(master=parent, text=text, command=command)
        else:
            b = Tk.Button(master=parent, text=text, padx=2, pady=2, command=command)
        b.pack(side=side)
        return b

    def _cancel_text(self):
        self._finish_textwindow()

    def _print_text(self):
        self.print_text()
        self._finish_textwindow()

    def load_textwindow(self,event):
        if event.canvas._tkcanvas != self.parent:
            raise RuntimeError, "Got invalid event!"
        
        self.event = event
        is_ax = (event.inaxes != None)
        (xpix, ypix) = self._disppix2screen(event.x, event.y)
        offset=5
        self.show_textwindow(xpix+offset,ypix+offset,enableaxes=is_ax)
        
    def show_textwindow(self,xpix,ypix,basetext=None,enableaxes=False):
        if not self.textwin: return
        self._reset_radio()
        if enableaxes: 
            self._enable_radio()
        self.textwin.deiconify()
        (w,h) = self.textwin.minsize()
        if w*h <= 1:
            self.textwin.minsize(width=self.textwin.winfo_width(),
                                 height=self.textwin.winfo_height())
            (w,h) = self.textwin.minsize()
        self.textwin.geometry("%sx%s+%s+%s"%(w,h,xpix,ypix))

    def _finish_textwindow(self):
        self._reset_radio()
        self._clear_textbox()
        self.textwin.withdraw()

    def _create_modmenu(self,master=None):
        if master:
            self.parent=master
        if not self.parent:
            return False
        menu = Tk.Menu(master=self.parent,tearoff=False)
        menu.add_command(label="Modify",command=self.modify)
        menu.add_command(label="Delete",command=self.delnote_dialog)
        return menu

    def show_modmenu(self,event):
        self.seltext = self._selected_text(event)
        if len(self.seltext) == 3:
            tkcanvas=event.canvas._tkcanvas
            xpos = tkcanvas.winfo_rootx() + int(event.x)
            ypos = tkcanvas.winfo_rooty() \
                   + tkcanvas.winfo_height() - int(event.y)
            self.menu.post(xpos,ypos)

    def modify(self):
        #print "Modify selected!!"
        textobj = self.seltext['textobj']
        (xtx, ytx) = textobj._get_xy_display()
        is_ax = (self.seltext['anchor'] != 'figure')
        if not is_ax:
            # previous anchor is figure
            pos = textobj.get_position()
            is_ax = (self._get_axes_from_pos(pos,self.canvas) != None)

        (xpix, ypix) = self._disppix2screen(xtx,ytx)
        offset = int(textobj.get_size())*2
        self.show_textwindow(xpix,ypix+offset,basetext=textobj.get_text(),\
                             enableaxes=is_ax)
        self._select_radio(self.seltext['anchor'])
        self._set_note(textobj.get_text())

    def delnote_dialog(self):
        remind = "Delete text?\n '"+self.seltext['textobj'].get_text()+"'"
        answer = tkMessageBox.askokcancel(parent=self.parent,title="Delete?",
                                          message=remind,
                                          default=tkMessageBox.CANCEL)
        if answer:
            self.delete_note()
        else:
            self.cancel_delete()


    def _disppix2screen(self,xpixd,ypixd):
        # calculate a pixel position form Upper-left of the SCREEN
        # from a pixel from Lower-left of the CANVAS (e.g., event.x/y)
        xpixs = self.parent.winfo_rootx() + xpixd
        ypixs = self.parent.winfo_rooty() + self.parent.winfo_height() \
               - ypixd
        return (int(xpixs), int(ypixs))
        
