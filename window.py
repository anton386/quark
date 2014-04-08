class Window(object):


    def __init__(self, allele_of_interest,
                 window_range,
                 increment_window):
        self.allele_of_interest = allele_of_interest
        if "-" in window_range:
            self.smallest_window = int(window_range.split("-")[0])
            self.largest_window = int(window_range.split("-")[1])
            self.window = None
        else:
            self.window = int(window_range)
        self.increment_window = increment_window
        self.windows = self._get_windows()
        self.num_of_windows = len(self.windows)

    def _get_windows(self):
        windows = []
        self.windows_dict = {}  # looks duplicated, but this is required for easier downstream retrieval

        if self.window:
            window_size = self.window
            windows.append((self.allele_of_interest - window_size/2,
                            self.allele_of_interest + window_size/2))
            self.windows_dict[window_size] = (self.allele_of_interest - window_size/2,
                                              self.allele_of_interest + window_size/2)
        else:
            window_size = self.smallest_window
            while window_size < self.largest_window:
                windows.append((self.allele_of_interest - window_size/2,
                                self.allele_of_interest + window_size/2))
                self.windows_dict[window_size] = (self.allele_of_interest - window_size/2,
                                                  self.allele_of_interest + window_size/2)
                window_size += self.increment_window

        return windows

    def get_window_start_by_size(self, size):
        return self.windows_dict[size][0]

    def get_window_end_by_size(self, size):
        return self.windows_dict[size][1]
