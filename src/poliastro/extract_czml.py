
class ExtractorCZML:
    """A class for extracting orbitary data to Cesium"""

    def __init__(self, orbit):
        """
        Orbital constructor

        :param orbit: poliastro.Orbit
            Orbit to be extracted
        """
        self.czml = dict()
        self.orbit = orbit

        self.start_epoch = ExtractorCZML.format_date(orbit.epoch.value)
        self.end_epoch = ExtractorCZML.format_date((orbit.epoch + orbit.period).value)

        self.init_czml()

        return None

    def parse_dict_tuples(self, path, tups):
        """
        Parameters
        ----------
        path : list (val)
            Dictionary path to insert to
        tups : list (val, val)
            Tuples to be assigned
        """

        # We only want to pass a reference czml, then by modifying our reference, we'll be modifying our base dictionary
        # which allows us to walk through a path of arbitrary length

        curr = self.czml
        for p in path:
            curr = curr.setdefault(p, {})
        for t in tups:
            curr[t[0]] = t[1]

    def init_czml(self):
        """
        Builds packet
        """
        self.parse_dict_tuples(["b0"], [("id", "document"), ("name", "simple"), ("version", "1.0")])
        self.parse_dict_tuples(["b0", "clock"], [("interval", self.start_epoch + '\\' + self.end_epoch),
                                                 ("currentTime", self.start_epoch), ("multiplier", 60),
                                                 ("range", "LOOP_STOP"),
                                                 ("step", "SYSTEM_CLOCK_MULTIPLIER")])

    def extract(self, ext_location):
        """
        Parameters
        ----------
        ext_location : str
            Path to extract your file to
        """

        return self.czml

    @staticmethod
    def format_date(date):
        """
        Parameters
        ----------
        date : str
            date of the form "yyyy-mm-dd hh:mm:ss.ssss"

        Returns
        -------
        formatted_date : str
            date of the form "yyyy-mm-ddThh:mm:ssZ"
        """
        return date[:10] + 'T' + date[11:-4] + 'Z'

