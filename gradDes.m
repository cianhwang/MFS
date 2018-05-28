function [A2 B2 E] = gradDes(A, B, Iv, In)
    global H_SOBEL
    H_SOBEL = fspecial('sobel');
    global H_LAPLACIAN
    H_LAPLACIAN = fspecial('laplacian', 0);
    global GRAD_IN
    GRAD_IN = gradient(In);
    global GRAD_IV
    GRAD_IV = gradient(Iv);
    THRESHOLD_E = 1e-1;
    THRESHOLD_L =1e-7;
    for i = 1:100
        [gA gB] = Eq11(A, B, Iv, In);
        LA = 0.5;
        E = energyCost(A, B, Iv, In);
        while(LA > THRESHOLD_L)
            dEA = energyCost(A-LA.*gA, B, Iv, In) - E;
            if (dEA > 0)
                LA = LA*0.7;
            else
                A = A - LA.*gA;
                break;
            end
            if (abs(dEA) < THRESHOLD_E)
                break;
            end
        end
        LB = 0.5;
        E = energyCost(A, B, Iv, In);
        while(LB > THRESHOLD_L)
            dEB = energyCost(A, B-LB.*gB, Iv, In) - E;
            if (dEB > 0)
                LB = LB*0.7;
            else
                B = B - LB.*gB;
                break;
            end
            if (abs(dEB) < THRESHOLD_E)
                break;
            end
        end
        if ((abs(dEA) < THRESHOLD_E) && (abs(dEB) < THRESHOLD_E))
            break;
        end
    end
   
    A2 = A;
    B2 = B;
    clear global
    
function [gA gB] = Eq11(A, B, Iv, In)
    global H_LAPLACIAN
    global GRAD_IN
    global GRAD_IV
    lap_A = imfilter(A, H_LAPLACIAN, 'replicate');
    lap_B = imfilter(B, H_LAPLACIAN, 'replicate');
    gA = 2*In.*(A.*In+B-Iv)+2*GRAD_IN.*...
        (A.*GRAD_IN-GRAD_IV)+2*0.1.*lap_A;
    gB = 2*(A.*In+B-Iv)+2*0.1*lap_B;

function E = energyCost(Av, Bv, Iv, In)
    global GRAD_IN
    global GRAD_IV
    grad_Av = gradient(Av);
    grad_Bv = gradient(Bv);
    E = sum(sum((Av.*In + Bv - Iv).^2 + ...
        (Av.*GRAD_IN - GRAD_IV).^2+...
        0.1 * (grad_Av.^2+grad_Bv.^2)));