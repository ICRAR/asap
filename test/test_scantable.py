import unittest
from asap import scantable, selector, rcParams
rcParams["verbose"] = False

class ScantableTest(unittest.TestCase):
    def setUp(self):

        self.st = scantable("data/MOPS.rpf", average=True)

    def test_init(self):
        st = scantable("data/MOPS.rpf", average=False)
        self.assertEqual(st.ncycle(), 32)
        st = scantable("data/MOPS.rpf", average=True)
        self.assertEqual(st.ncycle(), 2)
        st = scantable("data/MOPS.rpf", unit="Jy")
        self.assertEqual(st.get_fluxunit(), "Jy")
        st = scantable("data/MOPS.rpf", unit="K")
        self.assertEqual(st.get_fluxunit(), "K")
        self.assertRaises(RuntimeError, scantable, "data/MOPS.rpf", unit="junk")
        st = scantable(["data/MOPS.rpf","data/MOPS.rpf"], average=False)
        self.assertEqual(st.nscan(), 4)

    def test_copy(self):
        st = self.st.copy()
        self.assertNotEqual(id(st), id(self.st))
        
    def test_drop_scan(self):
        st = self.st.drop_scan([1])
        self.assertEqual(st.nscan(), 1)

    def test_get_scan(self):
        st = self.st.get_scan([1])
        self.assertEqual(st.nscan(), 1)
        st = self.st.get_scan("Orion_SiO_R")
        self.assertEqual(st.get_sourcename()[-1], "Orion_SiO_R")
        self.assertEqual(st.nscan(), 1)

    def test_get_spectrum(self):
        spec = self.st.get_spectrum(0)
        self.assertAlmostEqual(max(spec), 215.279830933)

    def test_get_mask(self):
        spec = self.st.get_mask(0)
        self.assertEqual(len(spec), 4096)
        
    def test_set_spectrum(self):
        spec = [ 1.0 for i in range(self.st.nchan()) ]
        self.st.set_spectrum(spec, 0)
        spec1 = self.st.get_spectrum(0)
        self.assertAlmostEqual(max(spec1), 1.0)

    def test_selection(self):
        sel = selector()
        sel.set_polarisations("YY")
        self.st.set_selection(sel)
        self.assertEqual(self.st.getpolnos(), (1,))
        sel1 = self.st.get_selection()
        self.assertEqual(sel1.get_pols(), [1])

    def test_stats(self):
        stats = { 'min': 113.767166138,
                  'max':215.279830933, 'sumsq':128759200.0,
                  'sum':720262.375, 'mean':175.845306396,
                  'var':513.95324707, 'stddev':22.6705360413,
                  'avdev':16.3966751099, 'rms':177.300170898, 
                  'median':182.891845703} 
        for k,v in stats.iteritems():
            sval = self.st.stats(stat=k)
            self.assertAlmostEqual(sval[0], v)
        msk = self.st.create_mask([0,100], [3900,4096])
        self.assertAlmostEqual(self.st.stats("sum", msk)[0], 35216.87890625)

    def test_get_column_names(self):
        cnames = ['SCANNO', 'CYCLENO', 'BEAMNO', 'IFNO', 
                  'POLNO', 'FREQ_ID', 'MOLECULE_ID', 'REFBEAMNO',
                  'TIME', 'INTERVAL', 'SRCNAME', 'SRCTYPE', 
                  'FIELDNAME', 'SPECTRA', 'FLAGTRA', 'TSYS',
                  'DIRECTION', 'AZIMUTH', 'ELEVATION', 
                  'PARANGLE', 'OPACITY', 'TCAL_ID', 'FIT_ID',
                  'FOCUS_ID', 'WEATHER_ID', 'SRCVELOCITY',
                  'SRCPROPERMOTION', 'SRCDIRECTION',
                  'SCANRATE']
        self.assertEqual(self.st.get_column_names(), cnames)


if __name__ == '__main__':
    unittest.main()

