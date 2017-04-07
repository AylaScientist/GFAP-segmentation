
function [writing] = Periodic_sampling_aut (imagefiles, nfiles, ind_path, project_path)


    %%%% Black rectangle %%%%
        a = 200; %High (100)
        b = 250; %Width (110)
        writing = 'black rectangle successfully created';
    %%%% Input of the image and number of subsets %%%%
    
    for j=1:nfiles
        image=imagefiles(j).name;
        imagen = imread([project_path,image]);
        writing = 'image files readed successfully';
        
        [A, B, c] = size (imagen); %A = height of the image, B = width of 
        % the image, c = colors


        %%%% Number of subsets %%%%
        n = ceil(A/(1.5*a)) * ceil(B/(1.5*b));
       
        %%%% RANDOM %%%%
        %integer random number between 0 and A or B 
        %(pseudorandom generation with normal distrbution)
        
        % x = floor (B*rand(1));
        % y = floor (A*rand(1));
        % Tx = B / b;
        % Ty = A / a;
        if A<B
            periodA = floor(A/sqrt(n));
            periodB = periodA;
        else
            periodB = floor(B/sqrt(n));
            periodA = periodB;
        end

        % AreaFraction = (a*b*n)/(A*B);


        Mask = zeros(A,B,c);
        x0 = floor(rand(1)*periodA/2);
        y0 = floor(rand(1)*periodB/2);

        xstart = x0;
        ystart = y0;
        while  (xstart < A-a) % proceed row-wise
            while (ystart < B-b)
                Mask(xstart+1:xstart+a, ystart+1:ystart+b,:) = 1;
                ystart = ystart+ periodB;
            end
            ystart = y0;    % when row ended, start again
            xstart = xstart + periodA;
        end

        MaskLabels = bwlabel(Mask(:,:,1));
        nlabels = max(MaskLabels(:));

        %figure(100), imshow(Mask.*imagen)
        figure(101), imshow(Mask)

        
        img_ab = zeros(a, b, 3, nlabels,'uint8');
        for i = 1:nlabels
            [rows, cols] = find(MaskLabels == i);
            img_ab(:,:,:,i)=imagen(rows(1):rows(end), cols(1):cols(end),:);
            %img_ab(:,:,:,i)=imagen(rows(1):rows(a), cols(1):cols(b),:);
            %clf(img_ab(:,:,:,i))
            imwrite(img_ab(:,:,:,i),[ind_path,'figure',(num2str(i*j)),'.jpg']);
        end
        message = ['SuperImage number ', (num2str(j)), 'finished']
    end
    clearvars img_ab;
    writing = 'Meander sampling: success'
end

