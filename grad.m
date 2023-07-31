function gradP = grad(ncellx,ncellr,dx,dr,p)
for i = 1:ncellr
    for j = 1:ncellx
        if j>1 && j< ncellx
            gradPx(i,j) = (p(i,j+1) - p(i,j-1))/2/dx;
        elseif j == 1
            gradPx(i,j) = (p(i,j+1) - p(i,j))/2/dx;
        elseif j == ncellx
            gradPx(i,j) = (p(i,j)-p(i,j-1))/2/dx;
        else
            error( 'Error occurred.\n not expected!');
        end

        if i>1 && i<ncellr
            gradPr(i,j) = (p(i+1,j) - p(i-1,j))/2/dr;
        elseif i == 1
            gradPr(i,j) = (p(i+1,j) - p(i,j))/2/dr;
        elseif i == ncellr
            gradPr(i,j) = (p(i,j) - p(i-1,j))/2/dr;
        else
            error( 'Error occurred.\n not expected!');
        end
        gradP(i,j,:) = [gradPx(i,j),gradPr(i,j)];
    end
end
end
