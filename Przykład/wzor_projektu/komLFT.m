function varargout = komLFT(varargin)
% KOMLFT MATLAB code for komLFT.fig
%Program symuluj�cy system komunikacyjny z modulacj� LFM 
%Przy uruchomieniu nale�y koniecznie wpisa� wszystkie nastawy i u�y� wszystkich element�w
%GUI, nawet gdy s� ustawione w�a�ciwie (SNR, pr�bkowanie)!
%Nast�pnie mo�na je zmienia� w dowolnej kolejno�ci.

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @komLFT_OpeningFcn, ...
                   'gui_OutputFcn',  @komLFT_OutputFcn, ...
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
% End initialization code - DO NOT EDIT

% --- Executes just before komLFT is made visible.
function komLFT_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = komLFT_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_Callback(hObject, eventdata, handles)
handles.f=str2double(get(hObject,'String')); 
guidata(hObject, handles);
function f_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function B_Callback(hObject, eventdata, handles)
handles.B=str2double(get(hObject,'String')) ;
guidata(hObject, handles);
function B_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function T_Callback(hObject, eventdata, handles)
handles.T=str2double(get(hObject,'String')); 
guidata(hObject, handles);
function T_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function SNR_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');
function SNR_Callback(hObject, eventdata, handles)
handles.SNR=SNR;
guidata(hObject, handles);

function sn_Callback(hObject, eventdata, handles);
n=get(hObject,'Value');
m=round((n-0.5)*40);
SN=num2str(m);
handles.snr=str2double(SN);
set(handles.SNR,'string',SN);
guidata(hObject, handles);
function sn_CreateFcn(hObject, eventdata, handles);

function x_Callback(hObject, eventdata, handles)
handles.x= get(hObject,'String'); 
guidata(hObject, handles);
function x_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function pz_Callback(hObject, eventdata, handles)
z=get(hObject,'value') ;
set(handles.pk,'value',0);
if z==1
handles.p=1;             %Pr�bkowanie zwyk�e
else
handles.p=0;            %Pr�bkowanie kwadraturowe
end
guidata(hObject, handles);

function pk_Callback(hObject, eventdata, handles)
z=get(hObject,'Value'); 
set(handles.pz,'value',0);
if z==1
handles.p=0;
else
handles.p=1;    
end
guidata(hObject, handles);

function y_Callback(hObject, eventdata, handles)
handles.y=y;
guidata(hObject, handles);
function y_CreateFcn(hObject, eventdata, handles)

function error_Callback(hObject, eventdata, handles)
handles.error=error;
guidata(hObject, handles);
function error_CreateFcn(hObject, eventdata, handles)

%Start do program oblLFM.m 
function pushbutton1_Callback(hObject, eventdata, handles)
%Nastawy
f=10^3*handles.f;       %Cz�stotliwo�� [Hz]
B=10^3*handles.B;   %Pasmo {Hz]
T=10^-3*handles.T;   %Czas trwania bitu [s]
snr=handles.snr;        %Stosunek sygna�u do szumu
x=handles.x;               %Tekst
p=handles.p;              %Wska�nik metody pr�bkowania

oblLFM
