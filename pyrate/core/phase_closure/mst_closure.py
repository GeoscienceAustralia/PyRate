# function[mloop] = mst_closure(ifglist,epochlist,iwrite)
# %====================================================================
# %function[mloop] = mst_closure(ifglist,epochlist,iwrite)
# %
# % find minimum independent loops using MST and Dijkstra algorithms
# %
# %INPUT:
# %  ifglist: interferogram list, see getnml.m for its structure
# %  epochlist: epoch list
# %  iwrite: output network (default 1: output closure; 0: not)
# %
# %OUTPUT:
# %  mloop: structure of the loops, including epochs, cost, and ifgs of each loop
# %
# % Hua Wang, 31/07/2011
# %
# % 10/08/2012 HW: fix a bug so that redundant loops will not be reused
# %====================================================================
# if nargin<3
#   iwrite=1;
# end

# %--------------------------------------------------------
# %1. calculate independent network number and MST
# %--------------------------------------------------------
# nimages = length(epochlist.date);
# nifgs=length(ifglist.nml);
#
# %initial connceted component
# connect = eye(nimages);
#
# %cost matrix
# costmat = zeros(nimages)+inf;

# %find mst
# mstlist=[];
# for i=1:nifgs
#   mas=ifglist.masnum(i);
#   slv=ifglist.slvnum(i);
#
#   %locate master/slave image
#   locmas = find(connect(:,mas)==1);
#   locslv = find(connect(:,slv)==1);
#
#   %set cost for each ifg in mstlist as 1 (nanfrac??)
#   costmat(mas,slv)=1;
#   costmat(slv,mas)=1;
#
#   if locmas~=locslv
#     mstlist=[mstlist;i];
#     connect(locmas,:)=connect(locmas,:)+connect(locslv,:); %add slave to MST
#     connect(locslv,:)=[];                                  %remove slave from independent component list
#   end
# end

def calc_independent_networ_number_and_mst():

    return


# %count isolated trees
# cnt=sum(connect');
# %remove missing images
# connect(cnt==1,:)=[];
# ntrees = size(connect,1);
# %redundant ifg list
# rdtlist = [1:nifgs];
# rdtlist(mstlist)=[];
# nrdt=length(rdtlist);
# fprintf('%d trees and %d redundant ifgs generated using MST algorithm\n',ntrees,nrdt);
#
# %--------------------------------------------------------------
# %2. Dijkstra algorithm for each redundant ifg
# %--------------------------------------------------------------
# for irdt=1:nrdt
#
#   mas=ifglist.masnum(rdtlist(irdt));
#   slv=ifglist.slvnum(rdtlist(irdt));
#
#   %find epoch number of current network
#   for i=1:ntrees
#     if connect(i,mas)>0 || connect(i,slv)>0
#       nimages_cur=nnz(connect(i,:));
#       break;
#     end
#   end
#
#   s(1:nimages) = 0;              %s, vector, set of visited vectors
#   dist(1:nimages) = inf;         % it stores the shortest distance between the source node and any other node;
#   prev(1:nimages) = nimages+1;   % Previous node, informs about the best previous node known to reach each network node
#   dist(mas) = 0;
#
#   %break direct connections of all following redundant ifgs, revised on 10/08/2012
#   icost=zeros(nrdt-irdt+1,1);
#   for iflordt=irdt:nrdt
#     mas_iflo=ifglist.masnum(rdtlist(iflordt));
#     slv_iflo=ifglist.slvnum(rdtlist(iflordt));
#     icost(iflordt)=costmat(mas_iflo,slv_iflo);
#     costmat(mas_iflo,slv_iflo)=inf;
#     costmat(slv_iflo,mas_iflo)=inf;
#   end
#
#   %using current image number here because it can't search through isolated networks
#   while sum(s)~=nimages_cur
#     %find candidates
#     candidate=[];
#     for i=1:nimages
#       if s(i)==0
#         candidate=[candidate dist(i)];
#       else
#         candidate=[candidate inf];
#       end
#     end
#
#     %new start from the minimum distance point
#     [u_index u]=min(candidate);
#     s(u)=1;
#
#     %calculate distance for new points connected with u
#     for i=1:nimages
#       if (dist(u)+costmat(u,i))<dist(i)
#         dist(i)=dist(u)+costmat(u,i);
#         prev(i)=u;
#       end
#     end
#   end
#
#   %back search for the epochs in the loop
#   sp = [slv];
#   while sp(1)~=mas
#     sp=[prev(sp(1)) sp];
#   end
#   mloop(irdt).epoch=[sp sp(1)];  %close the loop
#   mloop(irdt).cost = dist(slv);
#
#   %restore the original cost for this redundant ifg, revised on 10/08/2012
#   for iflordt=irdt:nrdt
#     mas_iflo=ifglist.masnum(rdtlist(iflordt));
#     slv_iflo=ifglist.slvnum(rdtlist(iflordt));
#     costmat(mas_iflo,slv_iflo)=icost(iflordt);
#     costmat(slv_iflo,mas_iflo)=icost(iflordt);
#   end
#
#   %determine sign of the ifgs
#   npts=length(sp);
#   for j=1:npts
#     for ii=1:nifgs
#       if (mloop(irdt).epoch(j)==ifglist.masnum(ii)) & (mloop(irdt).epoch(j+1)==ifglist.slvnum(ii))
#         mloop(irdt).ifg(j)=ii;
#         break;
#       elseif (mloop(irdt).epoch(j)==ifglist.slvnum(ii)) & (mloop(irdt).epoch(j+1)==ifglist.masnum(ii))
#         mloop(irdt).ifg(j)=-ii;
#         break;
#       end
#     end
#   end
# end
#
# %--------------------------------------------------------------
# %3. output the loops
# %--------------------------------------------------------------
# if iwrite==1
#   fid=fopen('network.list','wt');
#   for i=1:nrdt
#     npts=length(mloop(i).epoch);
#     for j=1:npts-1
#       fprintf(fid,'%s  ',char(ifglist.nml(abs(mloop(i).ifg(j)))));
#     end
#     fprintf(fid,'\n');
#   end
#   fclose(fid);
# end
