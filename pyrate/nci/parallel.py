from __future__ import print_function
"""
  Title: parallel.py
  
  Author:  Duncan Gray, Duncan.gray@ga.gov.au 

  Description: Parallel class to allow EQRM to run on a cluster.
  
  Version: $Revision: 1624 $  
  ModifiedBy: $Author: dgray $
  ModifiedDate: $Date: 2010-04-21 11:45:36 +1000 (Wed, 21 Apr 2010) $
  
  Copyright 2007 by Geoscience Australia
"""
import math

import socket
from numpy import arange


class Parallel(object):
    """ Parallelise EQRM so it can run on a cluster.

    Attributes:
    rank: What is the id of this node in the cluster.
    size: How many processors are there in the cluster.
    node: name of the cluster node.
    is_parallel: True if parallel is operational
    file_tag: A string that can be added to files to identify who wrote the
      file.      
    _make_block_file: Does this node have data to write to a block
      file?  WARNING This attribute is tighly coupled to calc_lo_hi.
      It is assuming calc_lo_hi is only called with one value.
      (Assumption is currently true)
      
    """
    def __init__(self, is_parallel=True):
        """
        Use is_parallel = False to stop parallelism, eg when running
        several scenarios.
        """
        
        if is_parallel is True:
            try:
                import pypar
            except ImportError:
                self._not_parallel()
            else:
                if pypar.size() >= 2:
                    self.rank = pypar.rank()
                    self.size = pypar.size()
                    self.node = pypar.get_processor_name()
                    self.is_parallel = True
                    self.file_tag = "pypar" + str(self.rank)
                    self.log_file_tag = "pypar" + str(self.rank)
                else:
                    self._not_parallel()
        else:
            self._not_parallel()
            
        # Some constants to identify messages
        self.load_event_set = 0
    
    def all_striped_indices(self, elements):
        """
        Return the indices for all nodes given the number of elements,
        using the striping pattern.
        e.g. 
        indices [array([ 0,  2,  4,  6,  8, 10, 12]), 
                 array([ 1,  3,  5,  7,  9, 11, 13])
        """
        all_elements = arange(elements)
        indices = []
        for node in range(self.size):
            indices.append(all_elements[node::self.size])
        return indices
        
    def calc_all_indices(self, elements):
        """
        A wrapper for the blocking algorithm used.
        """
        return self.all_striped_indices(elements)
        
    def striped_indices(self, elements):
        """
        Return the indices for the current node given the number of elements,
        using the striping pattern
        """
        return arange(elements)[self.rank::self.size]
    
    def calc_indices(self, elements):
        """
        A wrapper for the blocking algorithm used.
        """
        return self.striped_indices(elements)        
            
    def calc_lo_hi(self, elements):
        """
        Calculate the low index and the high index of length elements,
        so each node can work on a section of an array.

        Args:
          elements: Lenght of array/list etc.

        Return:
        
        """
        # floor - Returns the largest integral value
        # that is not greater than x.
        L = int(math.floor(1.0*elements/self.size))
        

        M = elements - self.size*L
        
        if (self.rank < M):
            lo = self.rank*L + self.rank
            hi = lo + L + 1
        else:
            lo = self.rank*L + M
            hi = lo + L

        self.lo = lo
        self.hi = hi
        if hi == lo:
            self._make_block_file = 0
        else:
            self._make_block_file = 1

        
    def _not_parallel(self):
        """
        Set the attributes if there is only one node.
        """
        self.rank = 0
        self.size = 1
        self.node = socket.gethostname() # The host name
        self.is_parallel = False
        self.file_tag = ''
        self.log_file_tag = '-0' # this is so there is always a log-0.txt file.
            
    def barrier(self):
        """
        Synchronisation point. Makes processors wait until all 
               processors have reached this point.
        """
        if self.is_parallel is True:
            import pypar
            pypar.barrier()
    
    def send(self, *args, **kwargs):
        """
        Wrapper for pypar.send
        """
        if self.is_parallel is True:
            import pypar

            pypar.send(*args, **kwargs)
            
    def receive(self, *args, **kwargs):
        """
        Wrapper for pypar.receive
        """
        if self.is_parallel is True:
            import pypar
            return pypar.receive(*args, **kwargs)
        else:
            return None
        
    def waitfor(self, msg, source):
        """
        Block on wait for the provided message from the given source
        """
        go = False if self.is_parallel else True
        while go is False:
            incoming = self.receive(source=source)
            if incoming == msg:
                go = True
            
    def notifyworkers(self, msg):
        """
        Send all nodes that aren't rank==0 msg 
        """
        if self.is_parallel is True:
            for node in range(self.size)[1:]:
                self.send(msg, node)
    
    def calc_num_blocks(self):
        """
        pre-req: calc_lo_hi has been calculated - and only calculated once!
        """
        if self.is_parallel is True:
            import pypar
            #print "synchronise self.rank", self.rank
            if self.rank == 0:
                calc_num_blocks = self._make_block_file
                for source in range(1, self.size):
                    #print "waiting.."
                    received = pypar.receive(source)
                    #print "received", received
                    calc_num_blocks += received
                return calc_num_blocks
            else:
                #print "sending from ", self.rank
                pypar.send(self._make_block_file, 0)
                #print "sent from ", self.rank
                
    def finalize(self):
        """
        End being parallel
        """
        if self.is_parallel is True:
            import pypar
            pypar.finalize()
