function symvfmm(k)
  
% Seems to break for k > 16 - Better way would be to break up by groups
  
  syms X Y Z
  % All 6 symmetry types
  sMat = [X Y Z; X Z Y; Y X Z; Y Z X; Z X Y; Z Y X];
  
  K = k*k*k;
  NK = k*(k+1)*(k+2)/6;		

  % Create arrays of symbolic monomials
  for s = 1:6
    % Create poly matrix
    bPolys = [1 X Y Z];
    for i=1:k-1
      for j = (1 + (i)*(i+1)*(i+2)/6):((i+1)*(i+2)*(i+3)/6 - (i+1))
        [i j j-(i*(i+1)/2)];
        bPolys(j) = sMat(s,1)*bPolys(j-(i*(i+1)/2));
      end
      for j = ((i+1)*(i+2)*(i+3)/6 - (i)):((i+1)*(i+2)*(i+3)/6 - 1)
        [i j j-(i*(i+3)/2)];
        bPolys(j) = sMat(s,2)*bPolys(j-(i*(i+3)/2));
      end
      idx = (i+1)*(i+2)*(i+3)/6;
      [idx i*(i+1)*(i+2)/6];
      bPolys(idx) = sMat(s,3)*bPolys((i*(i+1)*(i+2)/6));
    end
    if s == 1
      A = bPolys;
    else
      A = [A;bPolys];
    end
  end

  % First row contains proper order, so use it to find the permutation arrays
  for i=1:6
    for j = 1:NK
      S(i,j) = find(A(1,:) == A(i,j));
    end
  end

  % First row contains proper order, so use it to find the permutation arrays
  for i=1:6
    for j = 1:NK
      S(i,j) = find(A(1,:) == A(i,j));
    end
  end

  % offset for Matlab/C
  S = S-1;

  PN = 0;
  cnt = 0;
  for a = 1:-2:-1
    for b = 1:-2:-1
      for c = 1:-2:-1
        % Create poly matrix
        bPolys = [1];
        for i=1:k-1
          for j = (1 + (i)*(i+1)*(i+2)/6):((i+1)*(i+2)*(i+3)/6 - (i+1))
            [i j j-(i*(i+1)/2)];
            bPolys(j) = a*bPolys(j-(i*(i+1)/2));
          end
          for j = ((i+1)*(i+2)*(i+3)/6 - (i)):((i+1)*(i+2)*(i+3)/6 - 1)
            [i j j-(i*(i+3)/2)];
            bPolys(j) = b*bPolys(j-(i*(i+3)/2));
          end
          idx = (i+1)*(i+2)*(i+3)/6;
          [idx i*(i+1)*(i+2)/6];
          bPolys(idx) = c*bPolys((i*(i+1)*(i+2)/6));
        end
        if cnt == 0
          PN = bPolys;
        else
          PN = [PN;bPolys];
        end
        cnt = cnt+1;
      end
    end
  end

  
  
  
  fid = fopen('test.txt','wt');

  fprintf(fid, '#define P %d\n#define N %d\n\n', 1, -1);
  for i=1:NK
    fprintf(fid, '#define B%d %d\n', S(1,i), S(1,i));
  end
  fprintf(fid,'\n');
  
  fprintf(fid, 'int basisSyms[%d][%d] = {\n', 6, NK);
  for i = 1:6
    fprintf(fid, '{');
    for j = 1:NK
      fprintf(fid, 'B%d', S(i,j));
      if j ~= NK
        fprintf(fid,', ');
      end
    end
    fprintf(fid, '},\n');
  end
  fprintf(fid, '};\n\n\n');

  fprintf(fid, 'int symsPosNeg[%d][%d] = {\n', 6, NK);
  for i = 1:8
    fprintf(fid, '{');
    for j = 1:NK
      if (PN(i,j) == -1)
        fprintf(fid, 'N');
      else
        fprintf(fid, 'P');
      end
      if j ~= NK
        fprintf(fid,', ');
      end
    end
    fprintf(fid, '},\n');
  end
  fprintf(fid, '};\n');

  
  
  fclose(fid);

  