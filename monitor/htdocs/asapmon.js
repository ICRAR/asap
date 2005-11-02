function validateranges() {
    var fields = ["plotwindow","brangewindow"];
}

function baselineswitch(par) {
    var fields = ["polyorder","baselinerange"];
    for (i=0; i<fields.length; ++i) {
        var cont = document.getElementById(fields[i]);
        if (par.checked) {
            this.value = "True";
            cont.style.display = "block";
        } else {
            cont.style.display = "none";
            this.value = "False";
        }
    }
}

function clearunit() {
    var fields = ["plotwindow","brangewindow"];
    for (i=0; i<fields.length; ++i) {
        var cont = document.getElementById(fields[i]);
        cont.value = "";
    }
}

function unitswitch(unitval) {
    var fields = ["velframe","veldoppler","velrest"];
    for (i=0; i<fields.length; ++i) {
        var cont = document.getElementById(fields[i]);
        if (unitval == "channel") {
            cont.style.display = "none";
        } else {
            if (unitval == "GHz") {
             if (fields[i] == "velrest" || fields[i] == "veldoppler") cont.style.display = "none";
             else cont.style.display = "block";
            } else {
                cont.style.display = "block";
            }
        }
    }
    var ulblcont = document.getElementById("prangeunit");
    ulblcont.innerHTML = "";
    ulblcont.innerHTML = unitval;
    var ulblcont = document.getElementById("brangeunit");
    ulblcont.innerHTML = "";
    ulblcont.innerHTML = unitval;
    clearunit();
}


function listFiles(name, serviceURL) {
        var callback = (arguments.length == 3 ? arguments[2] : null);
        var request = new Ajax.Request(serviceURL, {
                parameters: "path=" + name,
                onComplete: function(transport) {
                        var resolutions = processResolutions(transport.responseXML.documentElement);
                        if (callback != null) {
                                callback.process(resolutions);
                        } else {
                                return resolutions;
                        }
                }
        });
}

function processResolutions(result) {
    //var results = result.getElementsByTagName("Listing");

    var fnames = result.getElementsByTagName("File");
    var files = [];
    if (fnames.length > 0 ) {
        for (i=0;i<fnames.length;++i) {
            files.push(fnames[i].firstChild.nodeValue)
        }
    }
    return files;
}


function insertFields() {
    var opts = document.getElementById("directory");
    var path = opts.selectedIndex;
    listFiles(path,"http://localhost/cgi-bin/asapmon/filelist.py", callbackHandler);
    return;
}


var callbackHandler = {
    process: function(parm) {
        var opts = document.getElementById("filelist");
        var fileopt = opts.options;
        //clear
        for (i=0;i<fileopt.length;++i) {
            fileopt[i] = null;
        }
        fileopt.length = 0;
        for (i=0;i<parm.length;++i) {
            var opt = document.createElement("option");
            opt.text = parm[i];
            opt.value = i;
            fileopt.add(opt);
        }
        // last file is selected
        if (fileopt.length > 0)
            opt.selected = fileopt.length-1;
    }
}
