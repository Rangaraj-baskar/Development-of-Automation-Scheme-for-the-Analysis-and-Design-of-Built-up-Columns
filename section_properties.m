%% Section properties 
% Moment of Inertia(I_xx, I_yy, I_xy), Area(area), Center of Gravity(C_G), Principal moment of Inertia (P_I_xx,
% P_I_yy,), Angle of inclination of principal axis with global cordinate
% axis(A_O_I) are found from Excel file given as Input
%
% *The input is given in the following format in the excel sheet:*
%  
%           columnns    n   x       y   t
%                       1	5       1	0.1
%                       2	5       0	0.1
%                       3	2.5     0	0.1
%                       4	0       0	0.1
%                       5	0       3	0.1
%                       6	0       6	0.1
%                       7	0       9	0.1
%                       8	2.5     9	0.1
%                       9	5       9	0.1
%                       10	5       8	
%
% Thickness is given for only 9 nodes as there are only (n-1) number of
% elements for every n number of nodes in the single section.\n
% The program finds the properties for both, the Back to back attached
% sections as well as box sections based on the input given for the left
% cross-section alone. 


%% The Main Function called 

function [I_xx, I_yy, I_xy, P_I_uu, P_I_vv, A_O_I, c_g_x, c_g_y, total_area, I_area, I_cg_x, I_cg_y, II_xx, II_yy, ...
            II_uu, II_vv, box_area, box_cg_x, box_cg_y, box_I_xx, box_I_yy, box_I_uu, box_I_vv, I_section_strength, box_section_strength ...
            ] = section_properties(filename)
    [~,~,data] = xlsread(filename);
    [s_L, x_coordinates, y_coordinates, row, column] = segment_length(data);
    slope = slope_calculator(x_coordinates, y_coordinates, row);
    [x_L_M_I, y_L_M_I, thickness] = local_M_I(s_L, slope, data, row, column);
    Area = area(s_L, thickness);
    total_area = sum(Area);
    [m_p_x, m_p_y] = mid_point(x_coordinates, y_coordinates, row);
    [c_g_x, c_g_y] = center_of_gravity(Area, m_p_x, m_p_y);
    [d_x, d_y] = distance_of_element_from_CG(m_p_x, m_p_y, c_g_x, c_g_y);
    [I_xx, I_yy, I_xy, P_I_uu, P_I_vv, A_O_I] = global_properties(x_L_M_I,y_L_M_I, Area, d_x, d_y);
    [I_area, I_cg_x, I_cg_y, II_xx, II_yy, II_uu, II_vv] = I_section_properties(total_area, x_coordinates, c_g_x, c_g_y, thickness, row, I_xx, I_yy);
    [box_area, box_cg_x, box_cg_y, box_I_xx, box_I_yy, box_I_uu, box_I_vv] = box_section_properties(total_area, x_coordinates, y_coordinates, c_g_x, c_g_y, I_xx, I_yy);
    I_section_strength = strength_of_I_section(total_area, I_xx, I_yy, I_area, II_xx, II_yy);
    box_section_strength = strength_of_box_section(total_area, I_xx, I_yy, box_area, box_I_xx, box_I_yy);
end

%% The Function to find out each segment length based on co-ordinates is given below

function [s_L, x_coordinates, y_coordinates, row, column] = segment_length(data)
    [row , column] = size(data);
    s_L = zeros(row-2,1);
    for column_index = 1 : column
        if strcmp(data(1, column_index),'x') || strcmp(data(1, column_index),'X')
            x_coordinates = double(string(data(2:end , column_index)));
        end
        if strcmp(data(1, column_index),'y') || strcmp(data(1, column_index),'Y') 
            y_coordinates = double(string(data(2:end , column_index)));
        end
    end
    for x = 1 : size(x_coordinates)-1
        s_L(x,1) = sqrt(abs((x_coordinates(x+1,1)- x_coordinates(x,1))^2 - (y_coordinates(x+1,1) - y_coordinates(x,1))^2));
    end
end
   
%% The Function to find slope of each segment.

function slope = slope_calculator(x_coordinates, y_coordinates, row)
    slope = zeros(row-1, 1);
    for x = 1 : size(x_coordinates) - 1
        slope(x,1) = sqrt(abs(x_coordinates(x+1, 1) - x_coordinates(x, 1)) .\ abs(y_coordinates(x+1, 1) - y_coordinates(x, 1)));
        slope(x,1) = atand(slope(x,1));
    end
end

%% The function for Calculation of Local Moment of Inertias of each segment.

function [x_L_M_I, y_L_M_I, thickness] = local_M_I(s_L, slope, data, row, column)
    x_L_M_I = zeros(row-2, 1);
    y_L_M_I = zeros(row-2, 1);
    for column_index = 1 : column
        if strcmp(data(1, column_index), 't') || strcmp(data(1, column_index), 'T')
            thickness = double(string(data(2:end , column_index)));
        end
    end
    for t = 1:(row-2)
        if slope(t,1) == 0 
            x_L_M_I(t,1) = (s_L(t,1)*(thickness(t,1).^3))./12;
            y_L_M_I(t,1) = ((s_L(t,1).^3)*thickness(t,1))./12;
        elseif slope(t,1) == 90 
            y_L_M_I(t,1) = (s_L(t,1)*(thickness(t,1).^3))./12;
            x_L_M_I(t,1) = ((s_L(t,1).^3)*thickness(t,1))./12;
        else
             I_xx = (s_L(t,1)*(thickness(t,1).^3))./12;
             I_yy = ((s_L(t,1).^3)*thickness(t,1))./12;
             x_L_M_I(t,1) = (I_xx + I_yy)/2 + ((I_xx - I_yy)/2)*cosd(2*slope(t,1));
             y_L_M_I(t,1) = (I_xx + I_yy)/2 + ((I_xx - I_yy)/2)*cosd(2*slope(t,1));
        end 
    end 
end

%% The function to calculate all the elements area 

function a = area(s_L, thickness)
    a = s_L.*thickness(1:size(s_L),1);
end

%% The function to calculate mid point of each individual segment

function [m_p_x, m_p_y] = mid_point(x_coordinates, y_coordinates, row)
    m_p_x = zeros(row-2, 1);
    m_p_y = zeros(row-2, 1);
    for row_index = 1 : (row - 2)
        m_p_x(row_index, 1) = (x_coordinates(row_index, 1) + x_coordinates(row_index + 1, 1))*0.5;
        m_p_y(row_index, 1 ) = (y_coordinates(row_index, 1) + y_coordinates(row_index + 1, 1))*0.5;
    end
end

%% The function to calculate the CG of complete segment

function [c_g_x, c_g_y] = center_of_gravity(A, m_p_x, m_p_y)
    c_g_x = sum(A.*m_p_x)/sum(A);
    c_g_y = sum(A.*m_p_y)/sum(A);
end

%% The function to calculate the distance along x axis and y axis of element from c_g

function [d_x, d_y] = distance_of_element_from_CG(m_p_x, m_p_y, c_g_x, c_g_y)
   d_x = ( m_p_x - c_g_x );
   d_y = ( m_p_y - c_g_y );
end

%% The function to calculate Global properties 

function [I_xx, I_yy, I_xy, P_I_uu, P_I_vv, A_O_I] = global_properties(x_L_M_I,y_L_M_I, A, d_x, d_y)
    I_yy = sum(x_L_M_I + A.*(d_x.^2));
    I_xx = sum(y_L_M_I + A.*(d_y.^2));
    I_xy = sum(A.*(d_x.*d_y));
    P_I_uu = (I_xx + I_yy).*0.5 + sqrt((((I_xx - I_yy)*0.5).^2) + I_xy.^2);
    P_I_vv = (I_xx + I_yy).*0.5 - sqrt((((I_xx - I_yy)*0.5).^2) + I_xy.^2);
    A_O_I = atand((-2*I_xy)/(I_xx - I_yy));
end

%% The function to calculate Back-Back attached section properties. 

function [I_area, I_cg_x, I_cg_y, II_xx, II_yy, II_uu, II_vv] =...
            I_section_properties(total_area, x_coordinates, c_g_x, c_g_y, thickness, row, I_xx, I_yy)
    I_area = 2*total_area;
    element = input("Enter the web element number for getting thickness of the Built up section web thickness : ");
    if element > (row-2)
        error("Enter a valid element number");
    else
        t = thickness(element,1);
    end
    I_cg_x = min(x_coordinates) - t/2;
    I_cg_y = c_g_y;
    II_xx = 2*I_xx;
    II_yy = 2*(I_yy + total_area*(I_cg_x - c_g_x)^2);
    II_uu = (II_xx + II_yy)*0.5 + sqrt((((II_xx - II_yy)*0.5)^2));
    II_vv = (II_xx + II_yy)*0.5 - sqrt((((II_xx - II_yy)*0.5)^2));
end

%% The Function to calculate Box properties.

function [box_area, box_cg_x, box_cg_y, box_I_xx, box_I_yy, box_I_uu, box_I_vv] = ...
            box_section_properties(total_area, x_coordinates, y_coordinates, c_g_x, c_g_y, I_xx, I_yy)
    box_area = 2*total_area;
    box_cg_x = (max(x_coordinates) + min(x_coordinates))*0.5;
    box_cg_y = (max(y_coordinates) + min(y_coordinates))*0.5;
    box_I_xx = 2*(I_xx + total_area*(box_cg_y - c_g_y)^2);
    box_I_yy = 2*(I_yy + total_area*(box_cg_x - c_g_x)^2);
    box_I_uu = (box_I_xx + box_I_yy)*0.5 + sqrt((((box_I_xx - box_I_yy)*0.5)^2));
    box_I_vv = (box_I_xx + box_I_yy)*0.5 - sqrt((((box_I_xx - box_I_yy)*0.5)^2));
end

%% The Function to find Strength of an I section 

function I_section_strength = strength_of_I_section(total_area, ...
            I_xx, I_yy, I_area, II_xx, II_yy)
    ri = sqrt(min(I_xx, I_yy)/total_area);
    ro = sqrt(min(II_xx, II_yy)/I_area);
    support = input("Enter the support type for I section('Fixed' or 'pinned') : ", 's');
    support = lower(support);
    if strcmp('fixed',support)==1
        k = 0.5;
    elseif strcmp('pinned',support)==1
        k = 1;
    else
        error("Enter either 'pinned' or 'fixed' properly");
    end
    L = input("Enter the length of the I section column : ");
    a = input("Enter the fastener spacing for I section column : ");
    modified_slenderness = sqrt((k*L/ro)^2 + (a/ri)^2);
    Fcre = (pi^2)/modified_slenderness;
    Fy = input("Enter the yield strength of steel used : ");
    lambda = sqrt(Fy/Fcre);
    if lambda <= 1.5
        Fn = ((0.658^(lambda^2))*Fy);
    else 
        Fn = ((0.877/(lambda^2))*Fy);
    end
    Pne = I_area*Fn;
    I_section_strength = Pne;
end

%% The function to calculate Strength of Box section

function box_section_strength = strength_of_box_section(total_area, ...
            I_xx, I_yy, box_area, box_I_xx, box_I_yy)
    E = 2.03e5;
    ri = sqrt(min(I_xx, I_yy)/total_area);
    ro = sqrt(min(box_I_xx, box_I_yy)/box_area);
    support = input("Enter the support type for box section ('Fixed' or 'pinned') : ", 's');
    support = lower(support);
    if strcmp('fixed',support)==1
        k = 0.5;
    elseif strcmp('pinned',support)==1
        k = 1;
    else
        error("Enter either 'pinned' or 'fixed' properly");
    end
    L = input("Enter the length of the box section column : ");
    a = input("Enter the fastener spacing for box section column : ");
    modified_slenderness = sqrt((k*L/ro)^2 + (a/ri)^2);
    Fcre = ((pi^2)*E)/modified_slenderness;
    Fy = input("Enter the yield strength of steel used : ");
    lambda = sqrt(Fy/Fcre);
    if lambda <= 1.5
        Fn = ((0.658^(lambda^2))*Fy);
    else 
        Fn = ((0.877/(lambda^2))*Fy);
    end
    Pne = box_area*Fn;
    box_section_strength = Pne;
end