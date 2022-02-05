# generalcase computes the solution h to the differetial equation Lh=0 mod m;
# flist is the list of fi, m is the degree of the expected solution h, hlo is the starting point of h;
generalcase:=proc(flist,x,m,hlo)
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

# trail will make multiple trails to calculate the solution h and calculate the average sparsity;
# m is the degree of the expected solution h, num is the number of total trails, tau is the sparsity bound, d is the exponents bound;
trial:=proc(m,num,l,tau,d)
   local i, h, flist,h0,sparslist,avg;
   sparslist:=[0$num];
   for i from 1 to num do
       flist:=structure_b_flist(l,x,tau,d);
       h0:=1; # starting point;
       h:=generalcase(flist,x,m,h0);
       sparslist[i]:=nops([coeffs(h,x)]);
   end do;
   #print(sparslist);
   avg:=convert(sparslist,`+`)/num;
   #printf("the avg sparsity of solution h is : %f\n",avg);
   return evalf(avg);  # evalf changes fraction to decimal;
end:

# plotl will calculate a list of avg sparsity based on a list of l values;
# m is the degree of the expected solution h, num is the number of total trials each l, tau is the sparsity bound, d is the exponents bound, len is the length of the llist;
plotl:=proc(m,num,tau,d,len)
   local llist, avglist,i;
   llist:=[0$len];
   avglist:=[0$len];
   # llist is assigned by i, which is different from structure_a;
   for i from 1 to len do
        llist[i]:=i;
        avglist[i]:=trial(m,num,i,tau,d);
   end do;
   print(llist);
   print(avglist);
   return llist,avglist;
end:

# plottau will calculate a list of avg sparsity based on a list of tau values;
# m is the degree of the expected solution h, num is the number of total trials each l, l is the order, d is the exponents bound, len is the length of the llist;
plottau:=proc(m,num,l,d,len)
   local taulist, avglist,i;
   taulist:=[0$len];
   avglist:=[0$len];
   # taulist is assigned by i, which is different from structure_a;
   for i from 1 to len do
        taulist[i]:=i;
        avglist[i]:=trial(m,num,l,taulist[i],d);
   end do;
   print(taulist);
   print(avglist);
   return taulist,avglist;
end:


