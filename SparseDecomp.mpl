# sparse_de_solver computes the solution h to the differetial equation Lh=0 mod m;
# flist is the list of fi, m is the degree of the expected solution h, hlo is the starting point of h;
sparse_de_solver:=proc(flist,x,m,hlo)
   local i,j,k,res,hk,l,h,htemp,hderi,temp;
   h:=hlo;
   l:= nops(flist)-1;
   # build up res;
   res := flist[1]*h;
   htemp :=h;
   for i from 2 to l+1 do
      hderi :=diff(htemp,x);
      res := res +flist[i]*hderi;
      htemp := hderi;
   end do;
   # main loop;
   while (ldegree(res)<m+1) do
      k :=ldegree(res)+l;
      hk := -((k-l)!*coeff(res,x,k-l))/(k!*coeff(flist[l+1],x,0));
      h := h+hk*x^k;
      for j from 1 to l+1 do
         res+= (k!/(k-j+1)!)*hk*flist[j]*x^(k-j+1);
      end do;
      res := expand(res);
   od;
   return h;
end:

# structure_a_flist generates a list of random polynomials fi;
# l is the order of the differential equation Lh=0, i.e. this function will generate f0,f1,...fl;
# tau is the sparsity bound of all fi, d is the exponents bound of all fi;
structure_a_flist:=proc(l,x,tau,d)
  local flist::list, i, supportlist, R, j,f0, taulist::list,Rcoeff,coefflist,term,expo;
  expo:=rand(1..d);
  R:=rand(1..tau-1);
  flist:=[0$l+1];
  taulist:=[0$l+1];
  # generate sparsity list
  for j from 1 to l+1 do
      taulist[j]:=R();
  end do;
  # generate fi list
  for i from 1 to l+1 do
      term:=taulist[i];
      flist[i]:=randpoly(x,degree=expo(),terms=term);
  end do;
  flist[l+1]:=expand(x*(flist[l+1])+R());  # make sure fl has nonzero constant;
  return flist;
end:

# structure_b_flist generates a list of fi such that f0=…x^(e1)+…+..x^(ek), f1=…x^(e1+1)+…+…x^(ek+1),  f2=…x^(e1+2)+…+..x^(ek+2)…;
# l is the order of the differential equation Lh=0, i.e. this function will generate f0,f1,...fl;
# tau is the sparsity bound of all fi, d is the exponents bound of all fi;
structure_b_flist:=proc(l,x,tau,d)
  local flist::list,explist::list,i,j,k,expo,coeff;
  expo:=rand(1..d); #d is the expo bound, set to equal to m;
  flist:=[0$l+1];
  # generate exponents of fl;
  for i from 1 to tau-1 do
     explist[i]:=expo();
  end do;
  for k from 1 to l+1 do
     for j from 1 to tau-1 do
        coeff:=rand(1..10);
        flist[k]+=coeff()*x^(explist[j]+k-1); # generate flist by shifting the exponents;
     end do;
  end do;
  flist[l+1]+=1; # make sure fl has nonzero constant;
  return flist;
end:

# structure_c_flist generates a list of fi such that fi = diff(f_(i+1),x);
# l is the order of the differential equation Lh=0, i.e. this function will generate f0,f1,...fl;
# tau is the sparsity bound of all fi, d is the exponents bound of all fi;
structure_c_flist:=proc(l,x,tau,d)
  local flist::list,explist::list,i,j,k,expo,coeff,ftemp;
  expo:=rand(1..d);
  flist:=[0$l+1];
  # generate a exponent list;
  for i from 1 to tau-1 do
     explist[i]:=expo();
  end do;
  # generate fl based on the exponent list;
  for j from 1 to tau-1 do
      coeff:=rand(1..10);
      flist[l+1]+=coeff()*x^(explist[j]+l);
  end do;
  # fi is the derivative of fi+1;
  for k from l to 1 by -1 do
      ftemp:=diff(flist[k+1],x);
      flist[k]:=ftemp;
  end do;
  flist[l+1]+=1; # make sure fl has nonzero constant;
  return flist;  
end:


# trail_a will make multiple trails for structure a to calculate the solution h and calculate the average sparsity;
# m is the degree of the expected solution h, num is the number of total trails, tau is the sparsity bound, d is the exponents bound;
trial_a:=proc(m,num,l,tau,d)
   local i, hlo, flist,h0,sparslist,avg;
   sparslist:=[0$num];
   for i from 1 to num do
       flist:=structure_a_flist(l,x,tau,d);
       h0:=1;  # assign starting point h0;
       hlo:=sparse_de_solver(flist,x,m,h0);
       sparslist[i]:=nops([coeffs(hlo,x)]); # get the sparsity of h;
   end do;
   avg:=convert(sparslist,`+`)/num; # calculate the avg sparsity;
   return evalf(avg); # evalf changes fraction to decimal;
end:

# trail_b will make multiple trails for structure b to calculate the solution h and calculate the average sparsity;
# m is the degree of the expected solution h, num is the number of total trails, tau is the sparsity bound, d is the exponents bound;
trial_b:=proc(m,num,l,tau,d)
   local i, h, flist,h0,sparslist,avg;
   sparslist:=[0$num];
   for i from 1 to num do
       flist:=structure_b_flist(l,x,tau,d);
       h0:=1; # starting point;
       h:=sparse_de_solver(flist,x,m,h0);
       sparslist[i]:=nops([coeffs(h,x)]);
   end do;
   #print(sparslist);
   avg:=convert(sparslist,`+`)/num;
   #printf("the avg sparsity of solution h is : %f\n",avg);
   return evalf(avg);  # evalf changes fraction to decimal;
end:

# trail_c will make multiple trails for structure c to calculate the solution h and calculate the average sparsity;
# m is the degree of the expected solution h, num is the number of total trails, tau is the sparsity bound, d is the exponents bound;
trial_c:=proc(m,num,l,tau,d)
   local i, h, flist,h0,sparslist,avg;
   sparslist:=[0$num];
   for i from 1 to num do
       flist:=structure_c_flist(l,x,tau,d);
       h0:=1; #starting point;
       h:=sparse_de_solver(flist,x,m,h0);
       sparslist[i]:=nops([coeffs(h,x)]);
   end do;
   #print(sparslist);
   avg:=convert(sparslist,`+`)/num;
   #printf("the avg sparsity of solution h is : %f\n",avg);
   return evalf(avg);  # evalf changes fraction to decimal;
end:

# plotl_a will calculate a list of avg sparsity based on a list of l values for structure a;
# m is the degree of the expected solution h, num is the number of total trials each l, tau is the sparsity bound, d is the exponents bound, len is the length of the llist;
plotl_a:=proc(m,num,tau,d,len)
   local llist, avglist,i;
   llist:=[0$len];
   avglist:=[0$len];
   # llist is assigned by i+1;
   for i from 1 to len do
        llist[i]:=i+1;
        avglist[i]:=trial_a(m,num,llist[i],tau,d);
   end do;
   print(llist);
   print(avglist);
   return llist,avglist;
end:

# plottau_a will calculate a list of avg sparsity based on a list of tau values fpr structure a;
# m is the degree of the expected solution h, num is the number of total trials each l, l is the order, d is the exponents bound, len is the length of the taulist;
plottau_a:=proc(m,num,l,d,len)
   local taulist, avglist,i;
   taulist:=[0$len];
   avglist:=[0$len];
   # taulist is assigned by i+1;
   for i from 1 to len do
        taulist[i]:=i+1;
        avglist[i]:=trial_a(m,num,l,taulist[i],d);
   end do;
   print(taulist);
   print(avglist);
   return taulist,avglist;
end:

# plotl_b will calculate a list of avg sparsity based on a list of l values for structure b;
# m is the degree of the expected solution h, num is the number of total trials each l, tau is the sparsity bound, d is the exponents bound, len is the length of the llist;
plotl_b:=proc(m,num,tau,d,len)
   local llist, avglist,i;
   llist:=[0$len];
   avglist:=[0$len];
   # llist is assigned by i, which is different from structure_a;
   for i from 1 to len do
        llist[i]:=i;
        avglist[i]:=trial_b(m,num,i,tau,d);
   end do;
   print(llist);
   print(avglist);
   return llist,avglist;
end:

# plottau_b will calculate a list of avg sparsity based on a list of tau values for structure b;
# m is the degree of the expected solution h, num is the number of total trials each l, l is the order, d is the exponents bound, len is the length of the taulist;
plottau_b:=proc(m,num,l,d,len)
   local taulist, avglist,i;
   taulist:=[0$len];
   avglist:=[0$len];
   # taulist is assigned by i, which is different from structure_a;
   for i from 1 to len do
        taulist[i]:=i;
        avglist[i]:=trial_b(m,num,l,taulist[i],d);
   end do;
   print(taulist);
   print(avglist);
   return taulist,avglist;
end:

# plotl_c will calculate a list of avg sparsity based on a list of l values for structure c;
# m is the degree of the expected solution h, num is the number of total trials each l, tau is the sparsity bound, d is the exponents bound, len is the length of the llist;
plotl_c:=proc(m,num,tau,d,len)
   local llist, avglist,i;
   llist:=[0$len];
   avglist:=[0$len];
   # llist is assigned by i, which is different from structure_a;
   for i from 1 to len do
        llist[i]:=i;
        avglist[i]:=trial_c(m,num,i,tau,d);
   end do;
   print(llist);
   print(avglist);
   return llist,avglist;
end:

# plottau_c will calculate a list of avg sparsity based on a list of tau values for structure c;
# m is the degree of the expected solution h, num is the number of total trials each l, l is the order, d is the exponents bound, len is the length of the taulist;
plottau_c:=proc(m,num,l,d,len)
   local taulist, avglist,i;
   taulist:=[0$len];
   avglist:=[0$len];
   # taulist is assigned by i, which is different from structure_a;
   for i from 1 to len do
        taulist[i]:=i;
        avglist[i]:=trial_c(m,num,l,taulist[i],d);
   end do;
   print(taulist);
   print(avglist);
   return taulist,avglist;
end:
