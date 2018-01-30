function varargout = FTSA(varargin)
%FTSA MATLAB code file for FTSA.fig
%      FTSA, by itself, creates a new FTSA or raises the existing
%      singleton*.
%
%      H = FTSA returns the handle to a new FTSA or the handle to
%      the existing singleton*.
%
%      FTSA('Property','Value',...) creates a new FTSA using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to FTSA_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      FTSA('CALLBACK') and FTSA('CALLBACK',hObject,...) call the
%      local function named CALLBACK in FTSA.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FTSA

% Last Modified by GUIDE v2.5 01-Nov-2017 14:40:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FTSA_OpeningFcn, ...
                   'gui_OutputFcn',  @FTSA_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FTSA is made visible.
function FTSA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for FTSA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FTSA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FTSA_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ExecuteMethod.
function ExecuteMethod_Callback(hObject, eventdata, handles)
% hObject    handle to ExecuteMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(handles.popupmenu1,'String'); 
popupmenu1value = contents{get(handles.popupmenu1,'Value')};
switch popupmenu1value
   case 'JPEG Compression with ZigZag'
       img = getimage(handles.axes1);
       if isempty(img)
           errordlg('Please upload an image first','Original Image is missing');
           return;
       end
       if size(img,1)~=size(img,2)
           img=imresize(img,[1024 1024]);
       end
       if size(img,3)==3
           img = rgb2gray(img);
       end
       dou_img = double(img);
       prompt = {'Enter Block Size','Enter Zigzag Coefficient'};
       dlg_title = 'JPEG Compression Variables';
       def = {'8','15'};
       temp = inputdlg(prompt,dlg_title,[1 50],def);
       if isempty(temp)
           % user pressed cancel
           return
       end
       [bs, col] = temp{:};
       if(~isempty(bs)&&~isempty(col))
           bs=str2num(bs);
           col=str2num(col);
           if rem(log2(bs),1)~=0
               errordlg('Please Enter Power of 2 for Block Size','Invalid Variable');
           elseif col<1
               errordlg('Please Enter ZigZag Coefficient for proper value','Invalid Variable');
           else
               fun = @(block_struct) dct2(block_struct.data);
               dct_img = blockproc(dou_img,[bs bs],fun);
               
               zzfun = @(block_struct) ZigZagscan(block_struct.data,col);
               zz_img = blockproc(dct_img,[bs bs],zzfun);
               
               ifun = @(block_struct) idct2(block_struct.data);
               idct_img = blockproc(zz_img,[bs bs],ifun);
               
               error_img = dou_img - idct_img;

               imshow(uint8(idct_img),'Parent',handles.axes2);
               title(handles.axes2,'Reconstruction Image');
               imshow(uint8(error_img),'Parent',handles.axes3);
               title(handles.axes3,'Error Image');
               
               MSE = immse(idct_img,double(img));
               PSNR = psnr(idct_img,double(img),255);
               
               fprintf('MSE = %0.2f\n',MSE); 
               fprintf('PSNR = %0.2f dB\n',PSNR);
           end
       else
           errordlg('Please complete all entries','Invalid Variable');
       end
   case 'JPEG Compression with Quality Factor' 
       img = getimage(handles.axes1);
       if isempty(img)
           errordlg('Please upload an image first','Original Image is missing');
           return;
       end
       if size(img,1)~=size(img,2)
           img=imresize(img,[1024 1024]);
       end
       if size(img,3)==3
           img = rgb2gray(img);
       end
       dou_img = double(img);
       prompt = {'Enter Block Size','Enter Quality'};
       dlg_title = 'JPEG Compression Variables';
       def = {'8','50'};
       temp = inputdlg(prompt,dlg_title,[1 50],def);
       if isempty(temp)
           % user pressed cancel
           return
       end
       [bs, q] = temp{:};
       if(~isempty(bs)&&~isempty(q))
           bs=str2num(bs);
           q=str2num(q);
           if rem(log2(bs),1)~=0
               errordlg('Please Enter Power of 2 for Block Size','Invalid Variable');
           elseif q<0 && q>100
               errordlg('Please Enter Quality from 1 to 100','Invalid Variable');
           else
               
               fun = @(block_struct) dct2(block_struct.data);
               img_DCT = blockproc(dou_img,[bs bs],fun);
               
               qfun = @(block_struct) round(block_struct.data/q);
               XQ = blockproc(img_DCT,[bs bs],qfun);

               [XQBitStream, ~] = ImageEncode(XQ);
               XQBitRate = length(XQBitStream)/prod(size(XQ));

               dqfun = @(block_struct) round(block_struct.data*q);
               XQ = blockproc(XQ,[bs bs],dqfun);

               ifun = @(block_struct) idct2(block_struct.data);
               idct_img = blockproc(XQ,[bs bs],ifun);

               imshow(uint8(idct_img),'Parent',handles.axes2);
               title(handles.axes2,'Reconstruction Image');
               fprintf('Bit rate of the compressed image = %0.2f bits/pixel\n',XQBitRate);
               
               error_img = dou_img - idct_img;
               
               imshow(uint8(error_img),'Parent',handles.axes3);
               title(handles.axes3,'Error Image');
               
               MSE = immse(idct_img,double(img));
               PSNR = psnr(idct_img,double(img),255);
               
               fprintf('MSE = %0.2f\n',MSE); 
               fprintf('PSNR = %0.2f dB\n',PSNR); 
           end
       else
           errordlg('Please complete all entries','Invalid Variable');
       end
   case 'Fast Fourier Transform'
       img = getimage(handles.axes1);
       if isempty(img)
           errordlg('Please upload an image first','Original Image is missing');
           return;
       end
       if size(img,3)==3
           img = rgb2gray(img);
       end
       douimg = im2double(img); %Transfer the image into double data type
       prompt = 'Enter a cut-off frequency for a low-pass filter:';
       dlg_title = 'Cut-off frequency';
       def = {'0.40'};
       answer = inputdlg(prompt,dlg_title,[1 50],def);
       if isempty(answer)
           % user pressed cancel
           return
       end
       cutoff = str2double(answer);
       % If loop that checks if K0 is correct (in range of frequencies)
       fftimg = fft2(douimg); % Take Fourier Transform 2D
       if (cutoff >= min(fftimg(:))) && (cutoff <= max(fftimg(:)))
           %% Processing
           % Fourier transform
           fftimg = fft2(douimg); % Take Fourier Transform 2D
           fftshiftimg = 20*log(abs(fftshift(fftimg))); % Shift center; get log magnitude
           % Application of low pass filter in reconstruction
           %Image dimensions 
           [row,col] = size(douimg); %[height, width]
           %Sampling intervals 
           dx = 1; 
           dy = 1; 
           %Characteristic wavelengths 
           KX = (mod(1/2 + (0:(col-1))/col, 1) - 1/2) * (2*pi/dx); 
           KY = (mod(1/2 + (0:(row-1))/row, 1) - 1/2) * (2*pi/dy); 
           [KX,KY] = meshgrid(KX,KY); 
           %Filter formulation 
           lpf = (KX.*KX + KY.*KY < cutoff^2); 
           %Filter Application 
           rec = ifft2(lpf.*fftimg);
           %% Results
           %Graphics
           %show FFT, magnitude
           mesh(fftshiftimg,'Parent',handles.axes2); 
           colormap(hot); % Plot Fourier Transform as function
           title(handles.axes2,'Hand-made low-pass filter');
           % Reconstruction after filtering in frequency domain
           imshow(rec,'Parent',handles.axes3);
           title(handles.axes3,'Filtered Image');
       else
           errordlg('Please enter valid entries and click ok','Invalid Variable');
       end
    case 'Add Noise by Filter'
        img = getimage(handles.axes1);
        if isempty(img)
            errordlg('Please upload an image first','Original Image is missing');
            return;
        end
        if size(img,3)==3
            img = rgb2gray(img);
        end
        prompt = 'Enter one of noise filters: Gaussian, Poisson, SaleAndPepper, Speckle';
        dlg_title = 'Add Noise';
        def = {'Gaussian'};
        answer = inputdlg(prompt,dlg_title,[1 50],def);
        if isempty(answer)
           % user pressed cancel
           return
        end
        if(~isempty(answer))
            if strcmp(answer, 'Gaussian')
                img_Gaussian = imnoise(img,'gaussian');
                imshow(img_Gaussian,'Parent',handles.axes1);
                title(handles.axes1,'Gaussian Image');
            elseif strcmp(answer, 'Poisson')
                img_Poisson = imnoise(img,'poisson');
                imshow(img_Poisson,'Parent',handles.axes1); 
                title(handles.axes1,'Poisson Image');
            elseif strcmp(answer, 'SaleAndPepper')
                img_SaleAndPepper = imnoise(img,'salt & pepper');
                imshow(img_SaleAndPepper,'Parent',handles.axes1); 
                title(handles.axes1,'salt & pepper Image');
            elseif strcmp(answer, 'Speckle')
                img_Speckle = imnoise(img,'speckle');
                imshow(img_Speckle,'Parent',handles.axes1);
                title(handles.axes1,'Speckle Image');
            else
                errordlg('Please enter valid noise filter','Invalid Variable');
            end
        else
            errordlg('Please enter noise filter','Empty Variable');
        end
    case 'Filter Denoise'
       img = getimage(handles.axes1);
       if isempty(img)
           errordlg('Please upload an image first','Original Image is missing');
           return;
       end
       if size(img,3)==3
           img = rgb2gray(img);
       end
       
       prompt = 'Enter one of denoise filters: low-pass, median, average, disk, gaussian, laplacian, motion, prewitt, sobel';
       dlg_title = 'Reduce Noise';
       def = {'gaussian'};
       answer1 = inputdlg(prompt,dlg_title,[1 50],def);
       
       if isempty(answer1)
           % user pressed cancel
           return
       end
       if(~isempty(answer1))
           if strcmp(answer1, 'low-pass')
               prompt = 'Enter a cut-off frequency for a low-pass filter:';
               dlg_title = 'Cut-off frequency';
               def = {'0.40'};
               cufoff = inputdlg(prompt,dlg_title,[1 50],def);
               cutoff = str2double(cufoff);
               % Generate the lowpass filter (order, cut-off frequency)
               lp = fir1(32,cutoff); 
               lp_2D = ftrans2(lp);  % Convert to 2-dimensions
               img_rep = imfilter(double(img),lp_2D,'replicate');
               img_lp = mat2gray(img_rep);
               imshow(img_lp,'Parent',handles.axes2); 
               title(handles.axes2,'Low Pass Noise Filter');
               mesh(img_lp,'Parent',handles.axes3);
               title(handles.axes3,'Frequency Spectrum');
           elseif strcmp(answer1, 'median')
               % Median Filter
               img_med = medfilt2(img);
               % Filtered image 1, lowpass filter
               imshow(img_med,'Parent',handles.axes2); 
               title(handles.axes2,'Median Noise Filter');
               mesh(img_med,'Parent',handles.axes3);
               title(handles.axes3,'Frequency Spectrum');
           elseif strcmp(answer1, 'average')
               prompt = 'Enter h of size hsize for a average filter:';
               dlg_title = 'Average Filter hsize';
               def = {'3'};
               hsize = inputdlg(prompt,dlg_title,[1 50],def);
               hsize = str2double(hsize);
               h = fspecial('average',hsize);
               img_avg = imfilter(img,h,'replicate');
               imshow(img_avg,'Parent',handles.axes2); 
               title(handles.axes2,'Average Filter');
               mesh(img_avg,'Parent',handles.axes3);
               title(handles.axes3,'Frequency Spectrum');
           elseif strcmp(answer1, 'disk')
               prompt = 'Enter a radius for a disk filter:';
               dlg_title = 'Disk Radius';
               def = {'5'};
               radius = inputdlg(prompt,dlg_title,[1 50],def);
               radius = str2double(radius);
               h = fspecial('disk',radius);
               img_disk = imfilter(img,h,'replicate');
               imshow(img_disk,'Parent',handles.axes2); 
               title(handles.axes2,'Disk Filter');
               mesh(img_disk,'Parent',handles.axes3);
               title(handles.axes3,'Frequency Spectrum');
           elseif strcmp(answer1, 'gaussian')
               prompt = {'Enter hsize for a gaussian filter:', 'Enter sigma for a gaussian filter:'};
               dlg_title = 'Gaussian Filter Variables';
               def = {'5','0.5'};
               temp = inputdlg(prompt,dlg_title,[1 50],def);
               [hsize,sigma] = temp{:};
               hsize = str2double(hsize);
               sigma = str2double(sigma);
               h = fspecial('gaussian',hsize,sigma);
               img_gaussian = imfilter(img,h,'replicate');
               imshow(img_gaussian,'Parent',handles.axes2); 
               title(handles.axes2,'Gaussian Filter');
               mesh(img_gaussian,'Parent',handles.axes3);
               title(handles.axes3,'Frequency Spectrum');
           elseif strcmp(answer1, 'laplacian')
               prompt = 'Enter alpha for a laplacian filter:';
               dlg_title = 'Laplacian Filter alpha';
               def = {'0.20'};
               alpha = inputdlg(prompt,dlg_title,[1 50],def);
               alpha = str2double(alpha);
               h = fspecial('laplacian',alpha);
               img_laplacian = imfilter(img,h,'replicate');
               imshow(img_laplacian,'Parent',handles.axes2); 
               title(handles.axes2,'Laplacian Filter');
               mesh(img_laplacian,'Parent',handles.axes3);
               title(handles.axes3,'Frequency Spectrum');
           elseif strcmp(answer1, 'log')
               prompt = {'Enter hsize for a log filter:', 'Enter sigmafor a log filter:'};
               dlg_title = 'Log Fliter Variables';
               def = {'5','0.5'};
               temp = inputdlg(prompt,dlg_title,[1 50],def);
               [hsize,sigma] = temp{:};
               hsize = str2double(hsize);
               sigma = str2double(sigma);
               h = fspecial('log',hsize,sigma);
               img_log = imfilter(img,h,'replicate');
               imshow(img_log,'Parent',handles.axes2); 
               title(handles.axes2,'Log Filter');
               mesh(img_log,'Parent',handles.axes3);
               title(handles.axes3,'Frequency Spectrum');
           elseif strcmp(answer1, 'motion')
               prompt = {'Enter len for a motion filter:', 'Enter theta for a motion filter:'};
               dlg_title = 'Motion Filter Variables';
               def = {'9','0'};
               temp = inputdlg(prompt,dlg_title,[1 50],def);
               [len,theta] = temp{:};
               len = str2double(len);
               theta = str2double(theta);
               h = fspecial('motion',len,theta);
               img_motion = imfilter(img,h,'replicate');
               imshow(img_motion,'Parent',handles.axes2); 
               title(handles.axes2,'Motion Filter');
               mesh(img_motion,'Parent',handles.axes3);
               title(handles.axes3,'Frequency Spectrum');
           elseif strcmp(answer1, 'prewitt')
               h = fspecial('prewitt');
               img_prewitt = imfilter(img,h,'replicate');
               imshow(img_prewitt,'Parent',handles.axes2); 
               title(handles.axes2,'Prewitt Filter');
               mesh(img_prewitt,'Parent',handles.axes3);
               title(handles.axes3,'Frequency Spectrum');
           elseif strcmp(answer1, 'sobel')
               h = fspecial('sobel');
               img_sobel = imfilter(img,h,'replicate');
               imshow(img_sobel,'Parent',handles.axes2); 
               title(handles.axes2,'Sobel Filter');
               mesh(img_sobel,'Parent',handles.axes3);
               title(handles.axes3,'Frequency Spectrum');
           else
               errordlg('Please enter valid denoise filter','Invalid Variable');
           end
       else
           errordlg('Please enter denoise filter','Empty Variable');
       end  
   otherwise
      fprintf('Please Select a method');
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in UploadImage.
function UploadImage_Callback(hObject, eventdata, handles)
% hObject    handle to UploadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,filepath]=uigetfile({'*.*','All Files'});
  img=imread([filepath,filename]);
  imshow(img,'Parent',handles.axes1);
  handles.img=img;
  guidata(hObject,handles);


% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes1,'reset');
cla(handles.axes2,'reset');
cla(handles.axes3,'reset');


% --- Executes on button press in Refresh.
function Refresh_Callback(hObject, eventdata, handles)
% hObject    handle to Refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = getimage(handles.axes1);
if isempty(temp)
    return;
end
img = handles.img;
imshow(img,'Parent',handles.axes1);
cla(handles.axes2,'reset');
cla(handles.axes3,'reset');
