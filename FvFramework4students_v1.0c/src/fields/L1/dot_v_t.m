function c = dot_v_t(a,na,b,nb,dim)
% a is a vector array of physical dimension dim
% na is the number of elements in na
% b is a tensor array of physical dimension dim
% nb is the number of elements in nb
% PRE: na == nb || na==1 || nb==1

   els = eldsize(1,dim);

   if na==1 && nb==1
      if dim == 1
         c = zeros(els,1);
         c(1) = a(1)*b(1);
         return
      elseif dim == 2
         c = zeros(els,1);
         c(1) = a(1)*b(1) + a(2)*b(2);
         c(2) = a(1)*b(3) + a(2)*b(4);
         return
      else
         c = zeros(els,1);
         c(1) = a(1)*b(1) + a(2)*b(2) + a(3)*b(3);
         c(2) = a(1)*b(4) + a(2)*b(5) + a(3)*b(6);
         c(3) = a(1)*b(7) + a(2)*b(8) + a(3)*b(9);
         return
      end   
   end
   
   
   
   if nb==1
      if dim == 1
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(1,ii)*b(1);
         end
         return
      elseif dim == 2
         c = zeros(els,na);
         for ii=1:na         
            c(1,ii) = a(1,ii)*b(1) + a(2,ii)*b(2);
            c(2,ii) = a(1,ii)*b(3) + a(2,ii)*b(4);
         end
         return
      else
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(1,ii)*b(1) + a(2,ii)*b(2) + a(3,ii)*b(3);
            c(2,ii) = a(1,ii)*b(4) + a(2,ii)*b(5) + a(3,ii)*b(6);
            c(3,ii) = a(1,ii)*b(7) + a(2,ii)*b(8) + a(3,ii)*b(9);
         end
         return
      end     
   end
   
   
   
   if na==1
      if dim == 1
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1)*b(1,ii);
         end
         return
      elseif dim == 2
         c = zeros(els,nb);
         for ii=1:nb         
            c(1,ii) = a(1)*b(1,ii) + a(2)*b(2,ii);
            c(2,ii) = a(1)*b(3,ii) + a(2)*b(4,ii);
         end
         return
      else
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1)*b(1,ii) + a(2)*b(2,ii) + a(3)*b(3,ii);
            c(2,ii) = a(1)*b(4,ii) + a(2)*b(5,ii) + a(3)*b(6,ii);
            c(3,ii) = a(1)*b(7,ii) + a(2)*b(8,ii) + a(3)*b(9,ii);
         end
         return
      end
   end
   
   
   
   if na==nb
      if dim == 1
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1,ii)*b(1,ii);
         end
         return
      elseif dim == 2
         c = zeros(els,nb);
         for ii=1:nb         
            c(1,ii) = a(1,ii)*b(1,ii) + a(2,ii)*b(2,ii);
            c(2,ii) = a(1,ii)*b(3,ii) + a(2,ii)*b(4,ii);
         end
         return
      else
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1,ii)*b(1,ii) + a(2,ii)*b(2,ii) + a(3,ii)*b(3,ii);
            c(2,ii) = a(1,ii)*b(4,ii) + a(2,ii)*b(5,ii) + a(3,ii)*b(6,ii);
            c(3,ii) = a(1,ii)*b(7,ii) + a(2,ii)*b(8,ii) + a(3,ii)*b(9,ii);
         end
         return
      end
   end
   
   
   c = [];
   

end