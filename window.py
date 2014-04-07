class Window(object):


    def __init__(self, allele_of_interest,
                 smallest_window, largest_window,
                 increment_window):
        self.allele_of_interest = allele_of_interest
        self.smallest_window = smallest_window
        self.largest_window = largest_window
        self.increment_window = increment_window
        self.windows = self._get_windows()
        self.num_of_windows = len(self.windows)

    def _get_windows(self):
        windows = []
        self.windows_dict = {}
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
