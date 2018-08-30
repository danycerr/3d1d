clear all
close all
fid = fopen ('../level2.frm');
txt = fgetl (fid);
dimensions = fscanf(fid, '%d');
nnu=dimensions(1)
nna=dimensions(2)
fclose (fid);
fid = fopen ('../level2.amg');
ia = fscanf(fid, '%d',nna+1);
ja = fscanf(fid, '%d',nnu);
a = fscanf(fid, '%f',nnu);
fclose (fid);
A=zeros(nna);
for i = 1:nna
	nonzero=ia(i+1)-ia(i);
	A(i,i) = a(ia(i));
	 for j  = 1:nonzero-1
	A(i, ja(ia(i) + j) ) = a(ia(i) +j); 		
	 end
end
spy(A)
