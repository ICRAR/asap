=====================
Scantable Description
=====================

Introduction
============

The Scantable is the ASAP single-dish data container. It was designed
based on the data abstraction of the PKSreader API.

Schema
======

The Scantable consist of a main (casacore) Table and several sub-tables which
are referenced via _id_.

----------
Main Table
----------

The main table consists of global data (keywords) and row-based data (columns)

Columns
-------

* **SPECTRA** - *Float(nchannel)*

    the spectral channel data.

* **FLAGTRA** - *uChar(nchannel)*

    the corresponding flags

* **TSYS** - *Float(nchannel)*

    the channel-based system temperature values

* **FLAGROW** - *uInt*

    spectrum based flags

* **TIME** - *MEpoch*

    The mjd time when the observation took place

* **INTERVAL** - *Double*

    the integration time in seconds

* **DIRECTION** - *MDirection*

    the associated direction on the sky

* **AZIMUTH, ELEVATION** - *Float*

    the azimuth/elevation when the spectrum was collected

* **SRCNAME** - *string*

    the name of the source observered

* **SRCTYPE** - *Int*

    the tyep of the source, i.e. indicating if it is an on source scan or
    off source. This will be used for calibration

* **SCANNO, BEAMNO, POLNO, IFNO, CYCLENO** - *uInt*

    These columns index (0-based) the respective values.

    * SCANNO: the number odf the scan. A scan is usually multiple integrations
              (cycles)

    * CYCLENO: the integration number within a scan (sub-scan?)

    * IFNO: the index of the IF (spectral window)

    * BEAMNO: the index of the beam (in a multibeam system)

    * POLNO: the index of the polarisation, e.g. XX=0, YY=1, Real(XY)=2,
             Imag(XY)=3

* **REFBEAMNO** - *Int*

    optional index of the reference beam in a multibeam obervation

* **FREQ_ID, MOLECULE_ID, TCAL_ID, FOCUS_ID, WEATHER_ID, FIT_ID**


----------
Sub-Tables
----------

FREQUENCY
---------

blah blah


========================
Mapping to other formats
========================

MS
==
