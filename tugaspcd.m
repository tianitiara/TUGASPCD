 function varargout = tugaspcd(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tugaspcd_OpeningFcn, ...
                   'gui_OutputFcn',  @tugaspcd_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code 


% --- Executes just before tugaspcd is made visible.
function tugaspcd_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = tugaspcd_OutputFcn(hObject, eventdata, handles) 



% --- Executes on button press in buka.
function buka_Callback(hObject, eventdata, handles)
global image;
gambar=uigetfile('*.jpg','Select Image file');
image=imread(gambar);
% handles.image=image;                  %untuk menyimpan nilai variabel
% guidata(hObject,handles);             %untuk simpan object
axes(handles.gambar1);                  %input nilai variabel pada axes1

imshow(image);                          %menampilkan gambar hasil browse

% --- Executes on button press in zoomout.
function zoomout_Callback(hObject, eventdata, handles)
global image;
%membuat matriks untuk gambar baru dengan baris dan kolom diperkecil sebanyak 2 kali lipat
newImage = zeros(round(size(image,1)/2), round(size(image,2)/2), 3);
m = 1; n = 1;
for i = 1:size(newImage,1)                  %perulangan baris
    for j = 1:size(newImage,2)              %perulangan kolom
        newImage(i,j,:) = image(m,n,:);     %mengeset nilai pixel pada gambar baru dengan pixel pada gambar
        n = round(n+2);                     %iterasi kolom
    end
    m = round(m+2);                         %iterasi baris
    n = 1;                                  %set ulang iterasi
end
newImage = uint8(newImage);                 %mengonvert matriks menjadi tipe data sebuah gambar
figure, imshow(newImage);                   %menampilkan hasil zoom out pada pop up 

% --- Executes on button press in zoomin.
function zoomin_Callback(hObject, eventdata, handles)
global image;
row = 2*size(image,1);        %menyiapkan baris gambar baru sebanyak 2 kali lipat dari gambar asli
column = 2*size(image,2);     %menyiapkan kolom gambar baru sebanyak 2 kali lipat dari gambar asli
%membuat matriks untuk gambar baru berdasarkan baris dan kolom yang telah diinisialisasi
newImage = zeros(row, column, 3); 
m = 1; n = 1;                 %iterasi baris & kolom untuk perulangan matriks baru
for i = 1:size(image,1)       % perulangan baris gambar
    for j = 1:size(image,2)   % perulangan kolom gambar
        % 4 baris kode dibawah ini adalah untuk mengeset nilai pixel pada 
        % pixel gambar dengan baris i  dan kolom j kepada 4 pixel yang sama
        % pada gambar baru
        newImage(m,n,:) = image(i,j,:);
        newImage(m,n+1,:) = image(i,j,:);
        newImage(m+1,n,:) = image(i,j,:);
        newImage(m+1,n+1,:) = image(i,j,:);
        n = n+2;   % iterasi tambah 2 karena setiap kolom pixel baru yg kosong berada di 2 angka setelah dari j
    end
    m = m+2;       % iterasi tambah 2 karena setiap baris pixel baru yg kosong berada di 2 angka setelah dari i
    n = 1;         % set ulang kolom kembali menjadi kolom pertama
end
newImage = uint8(newImage);  %konversi gambar baru menjadi tipe data uint8 (gambar)
figure, imshow(newImage);    %menampilkan hasil zoom in pada pop up 

% --- Executes on button press in Grayscale.
function Grayscale_Callback(hObject, eventdata, handles)
global image;
red = image(:,:,1);                           % memisahkan warna merah dari image
green = image(:,:,2);                         % memisahkan warna hijau dari image
blue = image(:,:,3);                          % memisahkan warna biru dari image
gray = 0.33*red+0.33*green+0.33*blue;         % merubah RGB ke Gray
axes(handles.gambar2);
imshow (gray)                                 % menampilkan hasil Grayscale


% --- Executes on button press in brightnestambah.
function brightnestambah_Callback(hObject, eventdata, handles)
global image;
global image2;
image2 = double(image) + 30;                    %menambahkan intensitas cahaya
image = double(image2);                         %memasukan variabel baru ke dalam variabel image
axes(handles.gambar2);                          %untuk menempatkan hasil image2 pada axes
imshow(uint8(image2));                          %menampilkan hasil Brightnes

% --- Executes on button press in brightneskurang.
function brightneskurang_Callback(hObject, eventdata, handles)
global image;
global image2;
image2 = double(image) - 30;                    %mengurangi intensitas cahaya
image = double(image2);                         %untuk menempatkan hasil image2 pada axes
axes(handles.gambar2);                          %memasukan variabel baru ke dalam variabel image
imshow(uint8(image2));                          %menampilkan hasil Brightnes

% --- Executes on button press in brightneskali.
function brightneskali_Callback(hObject, eventdata, handles)
global image;
global image2;
image2 = double(image) * 3;                     %mengkalikan intensitas cahaya
image = double(image2);                         %untuk menempatkan hasil image2 pada axes
axes(handles.gambar2);                          %memasukan variabel baru ke dalam variabel image
imshow(uint8(image2));                          %menampilkan hasil Brightnes

% --- Executes on button press in brightnesbagi.
function brightnesbagi_Callback(hObject, eventdata, handles)
global image;
global image2;
image2 = double(image) / 3;                     %menambahkan intensitas cahaya
image = double(image2);                         %megkalikan intensitas cahaya
axes(handles.gambar2);                          %memasukan variabel baru ke dalam variabel image
imshow(uint8(image2));                          %menampilkan hasil Brightnes


function verticalE_Callback(hObject, eventdata, handles)


function verticalS_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function verticalE_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function verticalS_Callback(hObject, eventdata, handles)

function verticalstart_DeleteFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function verticalstart_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function horizontalS_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function horizontalS_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function horizontalE_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function horizontalE_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in crop.
function crop_Callback(hObject, eventdata, handles)

global image;
global hasilcrop;
SV = str2num(get(handles.verticalS, 'String'));     %untuk menyimpan value kedalam variabel VS
SH = str2num(get(handles.horizontalS,'String'));    %untuk menyimpan value kedalam variabel HS
EV = str2num(get(handles.verticalE,'String'));      %untuk menyimpan value kedalam variabel VE
EH = str2num(get(handles.horizontalE, 'String'));   %untuk menyimpan value kedalam variabel VS

%inisialisasi gambar baru dengan batasan baris dan kolom yang telah ditentukan
hasilcrop = image(SV:EV, SH:EH);                    
hasilcrop(:,:,1) = image(SV:EV, SH:EH,1);           %menyimpan nilai pixel ke gambar baru
hasilcrop(:,:,2) = image(SV:EV, SH:EH,2); 
hasilcrop(:,:,3) = image(SV:EV, SH:EH,3); 
axes(handles.gambar2);                              %untuk menempatkan hasil pada axes
imshow(uint8(hasilcrop));                           %membuka gambar


% --- Executes on button press in histogram.
function histogram_Callback(hObject, eventdata, handles)

global image;
[baris, kolom, colormap] = size(image)
t1=1:256;                                           %variabel nilai pixel pada grafik R
n1=0:255;                                           %interval/range dari grafik R

t2=1:256;                                           %variabel nilai pixel pada grafik R
n2=0:255;                                           %interval/range dari grafik R

t3=1:256;                                           %variabel nilai pixel pada grafik R
n3=0:255;                                           %interval/range dari grafik R

count=0;                                            %variabel untuk menghitung nilai pixel kedalam grafik

%perulangan untuk grafik R
for z=1:256
    for i=1:baris
        for j=1:kolom
            if image(i,j,1)==z-1
                count=count+1;
            end
        end
    end
    t1(z)=count;
    count=0;
end

%perulangan untuk grafik G
for z=1:256
    for i=1:baris
        for j=1:kolom
            if image(i,j,2)==z-1
                count=count+1;
            end
        end
    end
    t2(z)=count;
    count=0;
end

%perulangan grafik B
for z=1:256
    for i=1:baris
        for j=1:kolom
            if image(i,j,3)==z-1
                count=count+1;
            end
        end
    end
    t3(z)=count;
    count=0;
end

figure, plot(n1,t1,'R',n2,t2,'G',n3,t3,'B');            %memunculkan grafik RGB 



% --- Executes on button press in histeq.
function histeq_Callback(hObject, eventdata, handles)

global image;
image = 0.2*image(:,:,1)+0.5*image(:,:,2)+0.3*image(:,:,3); %merubah GB ke grey

baris = size(image,1);                                      %set nilai baris kedalam variabel kolom
kolom= size(image,2);                                       %set nilai kolom kedalam variabel baris
jumlah=baris*kolom;

hasil=uint8(zeros(baris,kolom));                            %untuk gambar akhir setelah equalisasi
f=zeros(256,1);                                             % variabel awal dari nilai pixel+1
cdf=zeros(256,1);                                           %variabel yang nilainya bersifat kumulatif dari pdf
pdf=zeros(256,1);                                           %variabel yang menyimpan hasil dari variabel fungsi bagi jumlah
cum=zeros(256,1);                                           % variabel sementara yang bernilai kumulatif dari variabel fungsi 
out=zeros(256,1);                                           % variabel akhir yang akan disimpan kedalam variabel hasil

%perulangan untuk set nilai fungsi dan pdf
for i=1:baris
    for j=1:kolom
        nilaipixel=image(i,j);
        f(nilaipixel+1)=f(nilaipixel+1)/jumlah;
    end
end
sum=0;
L=255;

%perulangan untuk mendapatkan nilai pixel yang baru
for i=1:size(pdf)
    sum=sum+f(i);
    cum(i)=sum;
    cdf(i)=cum(i)/jumlah;
    out(i)=round(cdf(i)*L);
end

%perulangan untuk mengset variabel out kedalam variabel hasil
for i=1:baris
    for j=1:kolom
        hasil(i,j)=out(image(i,j)+1);
    end
end

t=1:256;
n=0:255;
count=0;

for z=1:256
    for i=1:baris
        for j=1:kolom
            if image(i,j)==z-1
                count=count+1;
            end
        end
    end
    t(z)=count;
    count=0;
end
imshow(hasil);
figure,plot(n,t);


% --- Executes on button press in blur.
function blur_Callback(hObject, eventdata, handles)

global image;                                           %pemanggilan variabel gambar pada axes / GUI
[m, n, o]=size(image);                                  %meng-get ukuran baris serta kolom dari gambar
newimage=double(image);                                 %inisialisasi matriks gambar untuk hasil yang baru
k = [1/8 1/8 1/8
    1/8 1/8 1/8
    1/8 1/8 1/8];                                       %inisialisasi matriks yang membuat gambar menjadi blur
for i=2 : m-2                                           %perulangan baris
    for j=2 : n-2                                       %perulangan kolom
        sum=newimage(i-1,j-1,:)*k(1,1,:)+newimage(i,j-1,:)*k(2,1,:)...
            + newimage(i+1,j-1,:)*k(3,1,:)+newimage(i-1,j,:)*k(1,2,:)...
            + newimage(i,j,:)*k(2,2,:)+newimage(i+1,j,:)*k(3,2,:)...
            +newimage(i-1,j+1,:)*k(1,3,:)+newimage(i,j+1,:)*k(2,3,:)...
            +newimage(i+1,j+1,:)*k(3,3,:);              %pengalian matriks k dengan gambar
        image(i-1,j-1,:)=sum;                           %hasil pixel yang baru tersimpan
    end
end

axes(handles.gambar2);                                   %untuk menempatkan hasil image pada axes
imshow(uint8(image));                                    %menampikan gambar hasil blur ke layar  


% --- Executes on button press in sharp.
function sharp_Callback(hObject, eventdata, handles)

global image;                                    %pemanggilan variabel gambar pada axes/GUI
[m, n, o]=size(image);                           %meng-get ukuran baris serta kolom dari gambar                
newimage=double(image);                          %inisialisasi matriks gambar untuk hasil yang baru
k = [0 -1 0
    -1 5 -1
    0 -1 0];                                     %inisialisasi matriks yang membuat gambar menjadi sharp/tajam
for i=2 : m-2                                    %perulangan baris
    for j=2 : n-2                                %perulangan kolom
        sum=newimage(i-1,j-1,:)*k(1,1,:)+newimage(i,j-1,:)*k(2,1,:)...
            + newimage(i+1,j-1,:)*k(3,1,:)+newimage(i-1,j,:)*k(1,2,:)...
            + newimage(i,j,:)*k(2,2,:)+newimage(i+1,j,:)*k(3,2,:)...
            +newimage(i-1,j+1,:)*k(1,3,:)+newimage(i,j+1,:)*k(2,3,:)...
            +newimage(i+1,j+1,:)*k(3,3,:);       %pengalian matriks k dengan gambar
        image(i-1,j-1,:)=sum;                    %hasil pixel yang baru tersimpan
    end
end

axes(handles.gambar2);

imshow(uint8(image));                           %menampilkan hasil sharp gambar ke layar

gambarsharp = uint8(image);
imwrite(gambarsharp,'hasil sharpen.jpg');


% --- Executes on button press in edge.
function edge_Callback(hObject, eventdata, handles)

global image;                   %pemanggilan variabel gambar pada axes/GUI
[m, n, o]=size(image);          %meng-get ukuran baris serta kolom dari gambar
newimage=double(image);         %inisialisasi matriks gambar untuk hasil yang baru
k = [1 1 1
    1 -13 1
    1 1 1];                     %inisialisasi matriks yang dapat mendektesi edge dari gambar
for i=2 : m-2                   %perulangan baris
    for j=2 : n-2               %perulangan kolom
        sum=newimage(i-1,j-1)*k(1,1)+newimage(i,j-1)*k(2,1)...
            + newimage(i+1,j-1)*k(3,1)+newimage(i-1,j)*k(1,2)...
            + newimage(i,j)*k(2,2)+newimage(i+1,j)*k(3,2)...
            +newimage(i-1,j+1)*k(1,3)+newimage(i,j+1)*k(2,3)...
            +newimage(i+1,j+1)*k(3,3);      %pengalian matriks k dengan gambar
        image(i-1,j-1)=sum;                 %hasil pixel deteksi edge tersimpan
    end
end

axes(handles.gambar2);
imshow(uint8(image));                       %menampilkan gambar dengan deteksi edge ke layar
gambaredge = uint8(image);
imwrite(gambaredge,'hasil edge.jpg');       %menyimpan gambar hasil deteksi edge


function Ra_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Ra_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Rb_Callback(hObject, eventdata, handles)

function Rb_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ga_Callback(hObject, eventdata, handles)

function Ga_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gb_Callback(hObject, eventdata, handles)

function Gb_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ba_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function Ba_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bb_Callback(hObject, eventdata, handles)
% hObject    handle to Bb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bb as text
%        str2double(get(hObject,'String')) returns contents of Bb as a double


% --- Executes during object creation, after setting all properties.
function Bb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in threshold.
function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
RAwal = str2double(get(handles.Ra, 'String'));
RAkhir = str2double(get(handles.Rb, 'String'));
GAwal = str2double(get(handles.Ga, 'String'));
GAkhir = str2double(get(handles.Gb, 'String'));
BAwal = str2double(get(handles.Ba, 'String'));
BAkhir = str2double(get(handles.Bb, 'String'));

[Tinggi, panjang, lebar]= size(image);
newimage= zeros(Tinggi,panjang,lebar);
for i = 1 : Tinggi
   for j=1 : panjang
        if image(i,j,1)>RAwal && image(i,j,1)<RAkhir && image(i,j,2)>GAwal && image(i,j,2)<GAkhir && image(i,j,3)>BAwal && image(i,j,3)<BAkhir
            newimage(i,j,:)=image(i,j,:);
        end
    end
end
axes(handles.gambar2);
b = uint8(newimage);
imshow(b);



function nilaiX_Callback(hObject, eventdata, handles)
% hObject    handle to nilaiX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nilaiX as text
%        str2double(get(hObject,'String')) returns contents of nilaiX as a double


% --- Executes during object creation, after setting all properties.
function nilaiX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nilaiX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nilaiY_Callback(hObject, eventdata, handles)
% hObject    handle to nilaiY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nilaiY as text
%        str2double(get(hObject,'String')) returns contents of nilaiY as a double


% --- Executes during object creation, after setting all properties.
function nilaiY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nilaiY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in seed.
function seed_Callback(hObject, eventdata, handles)
% hObject    handle to seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
masukan_valueX = str2double(get(handles.nilaiX,'String'));
masukan_valueY = str2double(get(handles.nilaiY,'String'));

nilai_threshold = 30;

[Tinggi,lebar, panjang]=size(image);
SegmentasiRegion = zeros(Tinggi, lebar, panjang);
matriksbaru = zeros(Tinggi,lebar, panjang);
seed= image(masukan_valueX,masukan_valueY);

while matriksbaru ~= 1
    for i = 1 : Tinggi
        for j = 1 : lebar
            if image(i,j)<nilai_threshold+seed && matriksbaru(i,j)~=1
                SegmentasiRegion(i,j)=nilai_threshold;
                matriksbaru(i,j)=1;
            end
        end
    end
    nilai_threshold = nilai_threshold + 30;
    for i = 1 : Tinggi
        for j = 1 : lebar
            if matriksbaru(i,j)~=1
                seed = image(i,j);
                break;
            end
        end
    end
end

axes(handles.gambar2);
b =uint8(SegmentasiRegion);
imshow(b);


% --- Executes on button press in erosi.
function erosi_Callback(hObject, eventdata, handles)

global image;
% image = handles.image;
image = im2bw(image);               %mengubah image dari grayscale ke dalam bentuk biner
newimage = imcomplement(image);     %proses negatif / complement pada gambar
[n,o]=size(newimage);               %get baris dan kolom gambar negatif
imageBerubah = double(newimage);    %inisialisasi variabel hasil dengan gambar, dengan tipe data double

for i=2 : n-2               %perulangan baris
    for j=2 : o-2           %perulangan kolom
        if newimage(i,j)==1 
            if imageBerubah(i,j-1)==0 || imageBerubah(i,j+1)==0 || imageBerubah(i-1,j)==0 || imageBerubah(i+1,j)==0
                imageBerubah(i,j)=2;
            end
        end
    end
end

imageresult = zeros(n,o);

for i=1 : n                %perulangan baris
    for j=1 : o            %perulangan kolom
        if imageBerubah(i,j)==1
            imageresult(i,j)=1;
        end
    end
end

axes(handles.gambar2);
imshow(imageresult);        %menampilkan hasil dari proses erosi ke layar



% --- Executes on button press in dilasi.
function dilasi_Callback(hObject, eventdata, handles)

global image;
image = im2bw(image);           %mengubah image dari grayscale ke dalam bentuk biner 
newimage = imcomplement(image); %proses negatif / complement gambar grayscale
[m, n]=size(newimage);          %get baris dan kolom pada gambar
imageBerubah = newimage;        %inisialisasi variabel hasilDilasi dengan gambar negatif

for i=2 : m-2                   %perulangan baris
    for j=2 : n-2               %perulangan kolom
        if newimage(i,j)==1     %jika nilai pixel = 1, 
            imageBerubah(i,j-1) = 1;
            imageBerubah(i,j)   = 1;
            imageBerubah(i,j+1) = 1;
            imageBerubah(i-1,j) = 1;
            imageBerubah(i+1,j) = 1;
        end
    end
end
gambargray = uint8(newimage);
axes(handles.gambar2);
imshow(imageBerubah);           %memunculkan hasil dari proses dilasi
imwrite(gambargray,'hasil dilasi.jpg'); %menyimpan hasil dari proses dilasi


% --- Executes on button press in image_compress.
function image_compress_Callback(hObject, eventdata, handles)

global image;
% image = handles.image;
imagelama = image;
[m,n,o]=size(imagelama);
array = double(imagelama);


for p = 1:3
    %ambil array dimensi ke-p
    arr = array(:,:,p);
    %buat array baru
    arraybaru = zeros(1,256);

    newimage = zeros(m,n);

    %melakukan kompresi dengan menghilangkan info yang penting
    for i= 1:m
        for j= 1:n
            if mod(double(arr(i,j)),2)==0
                newimage(i,j)=((double(arr(i,j))+1))/2;
            else
                newimage(i,j) = double(arr(i,j))/2;
            end
        end
    end
    
    gambarbaru(:,:,p)=newimage;
end

gambarbarulagi = uint8(gambarbaru);

figure
subplot(1,2,1),imshow(image);
axes(handles.gambar2);
imshow(gambarbarulagi);

imwrite(gambarbarulagi,'hasil kompresi.jpg');
