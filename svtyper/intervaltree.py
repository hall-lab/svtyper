from __future__ import print_function

# A very poorly implemented IntervalTree

# Based on reading:
# https://en.wikipedia.org/wiki/Interval_tree
# http://www.geeksforgeeks.org/interval-tree/
# http://bpbio.blogspot.com/2008/11/python-interval-tree.html

# A region is a tuple that consists of the following:
#    (sample_name, chromosome, position, left_position, right_position)

class IntervalTree(object):
    def __init__(self, regions):
        region_intervals = [ RegionInterval(r) for r in regions ]
        self.top_node = self.construct_tree(region_intervals)

    def construct_tree(self, intervals):
        if not intervals: return None

        center = self.get_center_interval(intervals)
        center_start = center.start

        center_intervals = []
        left_intervals = []
        right_intervals = []

        for k in intervals:
            if k.end < center_start:
                left_intervals.append(k)
            elif k.start > center_start:
                right_intervals.append(k)
            else:
                center_intervals.append(k)

        return Node(
            center_start,
            center_intervals,
            self.construct_tree(left_intervals),
            self.construct_tree(right_intervals)
        )

    def get_center_interval(self, intervals):
        sorted_intervals = sort_region_intervals(intervals)
        length = len(sorted_intervals)
        midpoint = int(length/2)
        mid_interval = sorted_intervals[midpoint]
        return mid_interval

    def search(self, begin, end):
        overlapping = self.top_node.find(begin, end)
        labels = [ o.label for o in overlapping ]
        return labels

class RegionInterval(object):
    def __init__(self, region):
        self.start = region[3]
        self.end = region[4]
        self.label = region

class Node:
    def __init__(self, center_start, center_intervals, left_node, right_node):
        self.center_start = center_start
        self.center_intervals = sort_region_intervals(center_intervals)
        self.left_node = left_node
        self.right_node = right_node

    def find(self, begin, end):
        overlapping = [ i for i in self.center_intervals
                          if i.end >= begin and i.start <= end ]

        if self.left_node and begin <= self.center_start:
            overlapping += self.left_node.find(begin, end)

        if self.right_node and end >= self.center_start:
            overlapping += self.right_node.find(begin, end)

        return overlapping

def sort_region_intervals(region_intervals):
    return sorted(region_intervals, key=lambda x: x.start)

if __name__ == "__main__":
    regions = [ ('sampleA', 'chr1', 50, 45, 55),
                ('sampleA', 'chr1', 60, 54, 80),
                ('sampleA', 'chr1', 75, 70, 80), ]

    tree = IntervalTree(regions)

    # find none
    test_interval = (20, 25)
    data = tree.search(*test_interval)
    print("{} -- {}".format(test_interval, data))

    # find one
    test_interval = (46, 50)
    data = tree.search(*test_interval)
    print("{} -- {}".format(test_interval, data))

    # find two
    test_interval = (51, 54)
    data = tree.search(*test_interval)
    print("{} -- {}".format(test_interval, data))

    # find three
    test_interval = (30, 90)
    data = tree.search(*test_interval)
    print("{} -- {}".format(test_interval, data))

