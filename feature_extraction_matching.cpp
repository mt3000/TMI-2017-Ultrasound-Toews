

/*
 
 Copyright (c) 2017, Matthew Toews
 All rights reserved.
 
 By downloading, copying, installing or using the software you agree to this license.
 If you do not agree to this license, do not download, install, copy or use the software.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 The views and conclusions contained in the software and documentation are those
 of the authors and should not be interpreted as representing official policies,
 either expressed or implied, of the FreeBSD Project.
 
 */



#include "opencv2/highgui/highgui.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/imgproc/imgproc.hpp"


#include <iostream>
#include <fstream>

#ifdef WIN32
  #include <time.h>
  #include <sys\timeb.h>
#endif


//
// normalize_img_grey()
//
// Normalize a greyscale image to unit length.
// Image could have an ROI
//
int
_normalize(
			  IplImage* img_in
			  )
{

	int iyOffset = 0;
	int ixOffset = 0;
	int iheight = img_in->height;
	int iweight = img_in->width;

	if( img_in->roi )
	{
		iyOffset = img_in->roi->yOffset;
		ixOffset = img_in->roi->xOffset;
		iheight = img_in->roi->height;
		iweight = img_in->roi->width;
	}

	double dSumSqr = 0;
	for( int r = 0; r < iheight; r++ )
	{
		for( int c = 0; c < iweight; c++ )
		{
			double dVal = (CV_IMAGE_ELEM(img_in,float,iyOffset+r,ixOffset+c));
			dSumSqr += dVal*dVal;
		}
	}

	if( dSumSqr > 0 )
	{
		double dNormFactor = 1.0/sqrt(dSumSqr);
		dSumSqr = 0;
		for( int r = 0; r < iheight; r++ )
		{
			for( int c = 0; c < iweight; c++ )
			{
				(CV_IMAGE_ELEM(img_in,float,iyOffset+r,ixOffset+c)) *= dNormFactor;
				double dVal = (CV_IMAGE_ELEM(img_in,float,iyOffset+r,ixOffset+c));

				dSumSqr += dVal*dVal;
			}
		}
	}
	else
	{
		// Set image to uniform

		double dValue = 1.0/sqrt( (double)(iheight*iweight) );

		for( int r = 0; r < iheight; r++ )
		{
			for( int c = 0; c < iweight; c++ )
			{
				(CV_IMAGE_ELEM(img_in,float,iyOffset+r,ixOffset+c)) = dValue;
				double dVal = (CV_IMAGE_ELEM(img_in,float,iyOffset+r,ixOffset+c));

				dSumSqr += dVal*dVal;
			}
		}
	}

	return 1;
}

class HoughCodeSimilarity
{
public:

	HoughCodeSimilarity(
											 )
	{
		m_fAngleToBooth = -1;
		m_fDistanceToBooth = -1;
		m_fAngleDiff = -1;
		m_fScaleDiff = -1;
	}

	~HoughCodeSimilarity(
											  )
	{
	}


	int
	SetHoughCode(
			const float &fPollRow,	// Poll booth row
			const float &fPollCol,	// Poll booth col
			const float &fPollOri,	// Poll booth orientation
			const float &fPollScl	// Poll booth scale
			)
	{
		// Set distance to poll booth
		float fDiffRow = fPollRow - m_Key_fRow;
		float fDiffCol = fPollCol - m_Key_fCol;
		float fDistSqr = fDiffCol*fDiffCol + fDiffRow*fDiffRow;
		m_fDistanceToBooth = (float)sqrt( fDistSqr );

		// Set angle between feature orientation and poll booth
		float fAngleToBooth = (float)atan2( fDiffRow, fDiffCol );
		m_fAngleToBooth = fAngleToBooth;

		// Set angular and orientation differences
		m_fAngleDiff = fPollOri - m_Key_fOri;
		m_fScaleDiff = fPollScl / m_Key_fScale;

		return 0;
	}

	int
	FindPollBooth(
							  const float &fRow,
							  const float &fCol,
							  const float &fOri,
							  const float &fScale,
							  float &fPollRow,		// Return row value
							  float &fPollCol,		// Return col value,
							  float &fPollOri,		// Return orientation value
							  float &fPollScl,		// Return scale value
							  int   bMirror	= 0// Match is a mirror of this feature
							  ) const
	{
		if( !bMirror )
		{
			float fOriDiff = fOri - m_Key_fOri;
			float fScaleDiff = fScale / m_Key_fScale;

			float fDist = m_fDistanceToBooth*fScaleDiff;
			float fAngle = m_fAngleToBooth + fOriDiff;

			fPollRow = (float)(sin( fAngle )*fDist + fRow);
			fPollCol = (float)(cos( fAngle )*fDist + fCol);

			fPollOri = fOri + m_fAngleDiff;
			fPollScl = m_fScaleDiff*fScale;
		}
		else
		{
			// Here, we consider a mirrored poll both
			float fOriDiff = fOri - (-m_Key_fOri);
			float fScaleDiff = fScale / m_Key_fScale;

			float fDist = m_fDistanceToBooth*fScaleDiff;
			float fAngle = (-m_fAngleToBooth) + fOriDiff;

			fPollRow = (float)(sin( fAngle )*fDist + fRow);
			fPollCol = (float)(cos( fAngle )*fDist + fCol);

			fPollOri = fOri - m_fAngleDiff;
			fPollScl = m_fScaleDiff*fScale;
		}

		return 0;
	}



	float m_fAngleToBooth;	// With just this, we can draw a line through to the booth
	float m_fDistanceToBooth;	// With just this, we can draw a circle through booth

	float m_fAngleDiff;		// Poll booth angle - key angle
	float m_fScaleDiff;		// Poll booth scale / key scale

	float m_Key_fRow;
	float m_Key_fCol;
	float m_Key_fOri;
	float m_Key_fScale;

};

class KeypointGeometry
{
public:
	KeypointGeometry(){};
	~KeypointGeometry(){};

	int iIndex0;
	int iIndex1;

	float m_Key_fRow;
	float m_Key_fCol;
	float m_Key_fOri;
	float m_Key_fScale;	
};


using namespace cv;
using namespace std;


#define PI 3.1415926535897
float
angle_radians(
			  float fAngle )
{
	return (2.0f*PI*fAngle) / 180.0f;
}

#define PI 3.141592653589793
#define SCALE_DIFF log((float)1.5f)
#define ORIENTATION_DIFF (20.0f*PI/180.0F)
//#define TRANSLATION_ERROR (30.0f/43.0f)
#define TRANSLATION_ERROR 4.00

int
compatible_poll_booths_line_segment(
					   float fTrainingCol,
					   float fTrainingRow,
					   float fTrainingOri,
					   float fTrainingScl,

					   float fPollCol,
					   float fPollRow,
					   float fPollOri,
					   float fPollScl,

						float fTranslationErrorThres = TRANSLATION_ERROR,
			float fScaleErrorThres = SCALE_DIFF,
			float fOrientationErrorThres = ORIENTATION_DIFF
				   )
{
	float fLogTrainingPollScl = log( fTrainingScl );
	float fLogPollScl = log( fPollScl );
	float fLogSclDiff = fabs( fLogTrainingPollScl - fLogPollScl );

	float fDiffCol = fPollCol - fTrainingCol;
	float fDiffRow = fPollRow - fTrainingRow;
	float fDistDiff = sqrt( fDiffCol*fDiffCol + fDiffRow*fDiffRow );

	float fOriDiff = fPollOri > fTrainingOri ?
		fPollOri - fTrainingOri : fTrainingOri - fPollOri;
	assert( fOriDiff >= 0 && fOriDiff < 2*PI );
	if( fOriDiff > (2*PI - fOriDiff) )
	{
		fOriDiff = (2*PI - fOriDiff);
	}
	assert( fOriDiff >= 0 && fOriDiff < 2*PI );

	// If the best match and this feature are within the match
	// hypothesis

	//if(		fDistDiff	< TRANSLATION_ERROR*fTrainingScl
	//	&&	fLogSclDiff	< SCALE_DIFF
	//	&&	fOriDiff	< ORIENTATION_DIFF
	//	)
	if(		fDistDiff	< fTranslationErrorThres*fTrainingScl
		&&	fLogSclDiff	< fScaleErrorThres
		&&	fOriDiff	< fOrientationErrorThres
		)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

//
// houghTransform()
//
// The hough transform here simply flags matches as null if they 
// do not agree with the global similarity transform.
//
int
houghTransform(
			   const vector<KeyPoint>& queryKeypoints,
                       const vector< vector<KeyPoint> > & trainKeypoints,
                       vector<DMatch> & matches,
					   int refx, int refy
					   )
{
	float fRefX;
	float fRefY;

	vector< KeypointGeometry > vecKG;
	//vecKG.resize( mmatches.size()*mmatches[0].size() );

	for( int i = 0; i < matches.size(); i++ )
	{
			if( matches[i].trainIdx < 0 || matches[i].imgIdx < 0 )
			{
				continue;
			}
			int iQIndex = matches[i].queryIdx;
			int iTIndex = matches[i].trainIdx;
			int iIIndex = matches[i].imgIdx;
			const KeyPoint &keyQ = queryKeypoints[iQIndex];
			const KeyPoint &keyT = trainKeypoints[iIIndex][iTIndex];

			HoughCodeSimilarity hc;

			hc.m_Key_fRow = keyQ.pt.y;
			hc.m_Key_fCol = keyQ.pt.x;
			hc.m_Key_fOri = angle_radians(keyQ.angle);
			hc.m_Key_fScale = keyQ.size;

			KeypointGeometry kg;

			hc.SetHoughCode( refy, refx, 0, 10 );
			hc.FindPollBooth( keyT.pt.y, keyT.pt.x, angle_radians(keyT.angle), keyT.size, kg.m_Key_fRow, kg.m_Key_fCol, kg.m_Key_fOri, kg.m_Key_fScale );
			kg.iIndex0 = i;
			vecKG.push_back( kg );
	}

	int iMaxIndex = -1;
	int iMaxCount = 0;
	for( int i = 0; i < vecKG.size(); i++ )
	{
		KeypointGeometry &kg1 = vecKG[i];
		int iCount = 0;
		for( int j = 0; j < vecKG.size(); j++ )
		{
			KeypointGeometry &kg2 = vecKG[j];
			if( compatible_poll_booths_line_segment(
					   kg1.m_Key_fCol,
					   kg1.m_Key_fRow,
					   angle_radians(kg1.m_Key_fOri),
					   kg1.m_Key_fScale,
					   kg2.m_Key_fCol,
					   kg2.m_Key_fRow,
					   angle_radians(kg2.m_Key_fOri),
					   kg2.m_Key_fScale ) )
			{
				iCount++;
			}
		}
		if( iCount > iMaxCount )
		{
			iMaxIndex = i;
			iMaxCount = iCount;
		}
	}

	for( int i = 0; i < matches.size(); i++ )
	{
		KeypointGeometry &kg1 = vecKG[iMaxIndex];
		for( int j = 0; j < vecKG.size(); j++ )
		{
			KeypointGeometry &kg2 = vecKG[j];
			if( compatible_poll_booths_line_segment(
					   kg1.m_Key_fCol,
					   kg1.m_Key_fRow,
					   angle_radians(kg1.m_Key_fOri),
					   kg1.m_Key_fScale,
					   kg2.m_Key_fCol,
					   kg2.m_Key_fRow,
					   angle_radians(kg2.m_Key_fOri),
					   kg2.m_Key_fScale ) )
			{
			}
			else
			{
				matches[j].imgIdx = -1;
				matches[j].queryIdx = -1;
				matches[j].trainIdx = -1;
			}
		}
	}

	return iMaxCount;
}

//
// houghTransform()
//
// Call with one vector of keypoints.
//
int
houghTransform(
			   const vector<KeyPoint>& queryKeypoints,
                       const vector<KeyPoint> & trainKeypoints,
                       vector<DMatch> & matches,
					   int refx, int refy
					   )
{
	vector< vector<KeyPoint> > vvtrainKeypoints;
	vvtrainKeypoints.push_back( trainKeypoints );
	return houghTransform( queryKeypoints, vvtrainKeypoints, matches, refx, refy );
}

//
// houghTransform()
//
// Call with a double vector of matches.
//
int
houghTransform(
			   const vector<KeyPoint>& queryKeypoints,
                       const vector< vector<KeyPoint> > & trainKeypoints,
                       vector< vector<DMatch> > & matches,
					   int refx, int refy
					   )
{
	vector<DMatch> vvMatches;

	for( int i = 0; i < matches.size(); i++ )
	{
		for( int j = 0; j < matches[i].size(); j++ )
		{
			if( matches[i][j].trainIdx >= 0 )
				vvMatches.push_back( matches[i][j] );
		}
	}

	return houghTransform( queryKeypoints, trainKeypoints, vvMatches, refx, refy );
}



const string defaultDetectorType = "SURF";
const string defaultDescriptorType = "SURF";
const string defaultMatcherType = "FlannBased";
const string defaultQueryImageName = "../../opencv/samples/cpp/matching_to_many_images/query.png";
const string defaultFileWithTrainImages = "../../opencv/samples/cpp/matching_to_many_images/train/trainImages.txt";
const string defaultDirToSaveResImages = "../../opencv/samples/cpp/matching_to_many_images/results";

void printPrompt( const string& applName )
{
	cout << "/*\n"
	<< " * This is a sample on matching descriptors detected on one image to descriptors detected in image set.\n"
	<< " * So we have one query image and several train images. For each keypoint descriptor of query image\n"
	<< " * the one nearest train descriptor is found the entire collection of train images. To visualize the result\n"
	<< " * of matching we save images, each of which combines query and train image with matches between them (if they exist).\n"
	<< " * Match is drawn as line between corresponding points. Count of all matches is equel to count of\n"
	<< " * query keypoints, so we have the same count of lines in all set of result images (but not for each result\n"
	<< " * (train) image).\n"
	<< " */\n" << endl;

    cout << endl << "Format:\n" << endl;
    cout << "./" << applName << " [detectorType] [descriptorType] [matcherType] [queryImage] [fileWithTrainImages] [dirToSaveResImages] [maskImage]" << endl;
    cout << endl;

    cout << "\nExample:" << endl
         << "./" << applName << " " << defaultDetectorType << " " << defaultDescriptorType << " " << defaultMatcherType << " "
         << defaultQueryImageName << " " << defaultFileWithTrainImages << " " << defaultDirToSaveResImages << "mask.pgm" << endl;
}

int
maskMatchesByTrainImgIdx( const vector<DMatch>& matches, int trainImgIdx, vector<char>& mask )
{
	int iCount = 0;
    mask.resize( matches.size() );
    fill( mask.begin(), mask.end(), 0 );
    for( size_t i = 0; i < matches.size(); i++ )
    {
        if( matches[i].imgIdx == trainImgIdx )
		{
            mask[i] = 1;
			iCount++;
		}
    }
	return iCount;
}

int
maskMatchesByMatchPresent( const vector<DMatch>& matches, vector<char>& mask )
{
	int iCount = 0;
    for( size_t i = 0; i < matches.size(); i++ )
    {
		if( matches[i].queryIdx <= 0 )
            mask[i] = 0;
		if( mask[i] > 0 )
		{
			iCount++;
		}
    }
	return iCount;
}

void readTrainFilenames( const string& filename, string& dirName, vector<string>& trainFilenames )
{
    const char dlmtr = '/';

    trainFilenames.clear();

    ifstream file( filename.c_str() );
    if ( !file.is_open() )
        return;

    size_t pos = filename.rfind(dlmtr);
    dirName = pos == string::npos ? "" : filename.substr(0, pos) + dlmtr;
    while( !file.eof() )
    {
        string str; getline( file, str );
        if( str.empty() ) break;
        trainFilenames.push_back(str);
    }
    file.close();
}

//
// readTrainTransforms_mha()
//
// Read in relevant lines from mha file, used with PLUS ultrasound data
//
void readTrainTransforms_mha(
							 const string& filename,
							 string& dirName,
							 vector<string>& trainFilenames,
							 vector <Mat>& trainImages,
							 vector< vector<float> >& trainImageData
							 )
{
#ifdef WIN32
    const char dlmtr = '\\';
#else
    const char dlmtr = '/';
#endif

    trainFilenames.clear();

    ifstream file( filename.c_str(), 0x21 );
    if ( !file.is_open() )
        return;

    size_t pos = filename.rfind(dlmtr);
    dirName = pos == string::npos ? "" : filename.substr(0, pos) + dlmtr;

	int iImgCols = -1;
	int iImgRows = -1;
	int iImgCount = -1;

	int iCurrImage = 0;

	// Vector for reading in transforms
	vector< float > vfTrans;
	vfTrans.resize(12);

	// Read until get dimensions
    while( !file.eof() )
    {
       string str; getline( file, str );
       if( str.empty() ) break;
	   char *pch = &(str[0]);
	   if( !pch )
		   return;

	   if( strstr( pch, "DimSize =" ) )
	   {
		   if( sscanf( pch, "DimSize = %d %d %d", &iImgCols, &iImgRows, &iImgCount ) != 3 )
		   {
			   printf( "Error: could not read dimensions\n" );
			   return;
		   }
	   }
	   if( strstr( pch, "ProbeToTrackerTransform =" )
		   || strstr( pch, "UltrasoundToTrackerTransform =" ) )
	   {
		   // Parse name and transform
		   // Seq_Frame0000_ProbeToTrackerTransform = -0.224009 -0.529064 0.818481 212.75 0.52031 0.6452 0.559459 -14.0417 -0.824074 0.551188 0.130746 -26.1193 0 0 0 1 

		   char *pcName = pch;
		   char *pcTrans = strstr( pch, "=" );
		   pcTrans[-1] = 0; // End file name string pcName
		   //pcTrans++; // Increment to just after equal sign

		   string filename = dirName + pcName + ".png";// + pcTrans;
		   trainFilenames.push_back( filename );

		   char *pch = pcTrans;

			for( int j =0; j < 12; j++ )
			{
				pch = strchr( pch + 1, ' ' );
				if( !pch )
					return;
				vfTrans[j] = atof( pch );
				pch++;
			}
			trainImageData.push_back( vfTrans );
	   }
	   if( strstr( pch, "ElementDataFile = LOCAL" ) )
	   {
		   // Done reading
		   break;
	   }
	}

	int iPosition = file.tellg();
	file.close();

	FILE *infile = fopen( filename.c_str(), "rb" );
	//fseek( infile, iPosition, SEEK_SET );
	char buffer[400];
	while( fgets( buffer, 400, infile ) )
	{
		if( strstr( buffer, "ElementDataFile = LOCAL" ) )
		{
			// Done reading
			break;
		}
	}
	int iPos2 = ftell( infile );

	unsigned char *pucImgData = new unsigned char[iImgRows*iImgCols];
	Mat mtImg(iImgRows, iImgCols, CV_8UC1, pucImgData);

	// Read & write images
	for( int i = 0; i < iImgCount; i++ )
	{
		Mat mtImgNew = mtImg.clone();
		fread( mtImgNew.data, 1, iImgRows*iImgCols, infile );
        trainImages.push_back( mtImgNew );
		//imwrite( trainFilenames[i], trainImages[i] );
		//imwrite( trainFilenames[i], mtImgNew );
    }

	for( int i = 0; i < iImgCount; i++ )
	{
		//imwrite( trainFilenames[i], trainImages[i] );
    }


	delete [] pucImgData;

	fclose( infile );

	printf( "Read images: %d\n", trainImages.size() );
	printf( "Read data  : %d\n", trainImageData.size() );
	printf( "Read names : %d\n", trainFilenames.size() );
    
}

//
// splitFile_mha()
//
// Split big mha file into several parts
//
void splitFile_mha(
							 const string& filename,
							 int iParts = 2
							 )
{
#ifdef WIN32
    const char dlmtr = '\\';
#else
    const char dlmtr = '/';
#endif

    ifstream file( filename.c_str(), 0x21 );
    if ( !file.is_open() )
        return;

	string dirName;

    size_t pos = filename.rfind(dlmtr);
    dirName = pos == string::npos ? "" : filename.substr(0, pos) + dlmtr;

	int iImgCols = -1;
	int iImgRows = -1;
	int iImgCount = -1;

	int iCurrImage = 0;

	// Vector for reading in transforms
	vector< float > vfTrans;
	vfTrans.resize(12);

	string fname1 = filename + "1.mha";
	string fname2 = filename + "2.mha";
	FILE *file1 = fopen( &fname1[0], "wb" );
	FILE *file2 = fopen( &fname2[0], "wb" );
	int iCount1;
	int iCount2;
	int bDone = 0;

	// For December AMIGO data (neuro-case/anonymized)
	// No intra-op, used for calibration
	int iStartRow = 124;
	int iStartCol = 265;
	int iRows = 640;
	int iCols = 690;

	// For January AMIGO data (20130108-us-anon\TrackedImageSequence_20130108_121528.mha)
	// 
	// For some reason there is an anoying shift ... 
	iStartRow = 140;
	iStartCol = 520;
	iRows = 750;
	iCols = 820;

	// For January AMIGO data (20130108-us-anon\TrackedImageSequence_20130108_132324.mha)
	// 
	// For some reason there is an anoying shift ... 
	iStartRow = 140;
	iStartCol = 520;

	// Andriy's prostate data december 2013 
	iStartRow = 0;
	iStartCol = 0;
	iRows = 700;
	iCols = 700;

	// Read until get dimensions
    while( !file.eof() && !bDone )
    {
       string str;
	   getline( file, str );
       if( str.empty() ) break;
	   char *pch = &(str[0]);
	   if( !pch )
		   return;

	   int iCurrFrame = -1;

	  if( iImgCount > -1 )
	  {
		  if( sscanf( pch, "Seq_Frame%d", &iCurrFrame ) != 1 )
		  {
			  iCurrFrame = -1;
		  }
		  else
		  {
			  printf( "Frame: %d\n", iCurrFrame );
		  }
	  }

	   if( strstr( pch, "DimSize =" ) )
	   {
		   if( sscanf( pch, "DimSize = %d %d %d", &iImgCols, &iImgRows, &iImgCount ) != 3 )
		   {
			   printf( "Error: could not read dimensions\n" );
			   return;
		   }
		   iCount1 = iImgCount / iParts;
		   iCount2 = iImgCount - iCount1;

		   // Special write operation
		   // Consider cropped images
		   fprintf( file1, "DimSize = %d %d %d",iCols,iRows,iCount1);
		   fprintf( file2, "DimSize = %d %d %d",iCols,iRows,iCount2);
		   continue;
	   }
	   //if( strstr( pch, "ProbeToTrackerTransform =" )
		  // || strstr( pch, "UltrasoundToTrackerTransform =" ) )
	   //{
		  // // Parse name and transform
		  // // Seq_Frame0000_ProbeToTrackerTransform = -0.224009 -0.529064 0.818481 212.75 0.52031 0.6452 0.559459 -14.0417 -0.824074 0.551188 0.130746 -26.1193 0 0 0 1 
	   //}
	   if( strstr( pch, "ElementDataFile = LOCAL" ) )
	   {
		   // Done reading
		   bDone = 1;
	   }

	   // Output
	   if( iCurrFrame < 0 )
	   {
		   fwrite( &str[0], 1, str.size(), file1 );fprintf( file1, "\n" );
		   fwrite( &str[0], 1, str.size(), file2 );fprintf( file2, "\n" );
	   }
	   else if( iCurrFrame < iCount1 )
	   {
		   fwrite( &str[0], 1, str.size(), file1 );fprintf( file1, "\n" );
	   }
	   else
	   {
		   fwrite( &str[0], 1, str.size(), file2 );fprintf( file2, "\n" );
	   }
	}

	int iPosition = file.tellg();
	file.close();

	//fclose( file1 );
	//fclose( file2 );

	FILE *infile = fopen( filename.c_str(), "rb" );
	fseek( infile, iPosition, SEEK_SET );

	unsigned char *pucImgData = new unsigned char[iImgRows*iImgCols];
	Mat mtImg(iImgRows, iImgCols, CV_8UC1, pucImgData);

	unsigned char *pucImgData2 = new unsigned char[iRows*iCols];
	Mat mtImg2(iRows, iCols, CV_8UC1, pucImgData2);

	// Read & write images
	for( int i = 0; i < iImgCount; i++ )
	{
		fread( mtImg.data, 1, iImgRows*iImgCols, infile );
		//imwrite( "image_big.png", mtImg );

		for( int r = 0; r < iRows; r++ )
		{
			memcpy( mtImg2.data+r*iCols, mtImg.data + (iStartRow + r)*iImgCols + iStartCol, iCols );
		}
		//imwrite( "image_small.png", mtImg2 );

		if( i < iCount1 )
		{
			for( int r = 0; r < iRows; r++ )
			{
				//fwrite( mtImg.data + (iStartRow + r)*iImgCols + iStartCol, 1, iCols, file1 );
				fwrite( mtImg2.data+r*iCols, 1, iCols, file1 );
			}
		}
		else
		{
			for( int r = 0; r < iRows; r++ )
			{
				//fwrite( mtImg.data + (iStartRow + r)*iImgCols + iStartCol, 1, iCols, file2 );
				fwrite( mtImg2.data+r*iCols, 1, iCols, file2 );
			}
		}
		//imwrite( trainFilenames[i], mtImg );
    }

	fclose( file1 );
	fclose( file2 );

	delete [] pucImgData;

	fclose( infile );
    
}

//
// convertToPng_mha()
//
// Convert mha file into png images
//
void
convertToPng_mha(
							 const string& filename,
							 const string& outputCode
							 )
{
#ifdef WIN32
    const char dlmtr = '\\';
#else
    const char dlmtr = '/';
#endif

	vector<string> trainFilenames;
	vector< vector<float> > trainImageData;

    trainFilenames.clear();

    ifstream file( filename.c_str() );
    if ( !file.is_open() )
        return;

    size_t pos = filename.rfind(dlmtr);
    string dirName = pos == string::npos ? "" : filename.substr(0, pos) + dlmtr;

	int iImgCols = -1;
	int iImgRows = -1;
	int iImgCount = -1;

	int iCurrImage = 0;

	// Vector for reading in transforms
	vector< float > vfTrans;
	vfTrans.resize(12);

	// Read until get dimensions
    while( !file.eof() )
    {
       string str; getline( file, str );
       if( str.empty() ) break;
	   char *pch = &(str[0]);
	   if( !pch )
		   return;

	   if( strstr( pch, "DimSize =" ) )
	   {
		   if( sscanf( pch, "DimSize = %d %d %d", &iImgCols, &iImgRows, &iImgCount ) != 3 )
		   {
			   printf( "Error: could not read dimensions\n" );
			   return;
		   }
	   }
	   if( strstr( pch, "ProbeToTrackerTransform =" )
		   || strstr( pch, "UltrasoundToTrackerTransform =" )
		   || strstr( pch, "ImageToCroppedImageTransform =" ) )
	   {
		   // Parse name and transform
		   // Seq_Frame0000_ProbeToTrackerTransform = -0.224009 -0.529064 0.818481 212.75 0.52031 0.6452 0.559459 -14.0417 -0.824074 0.551188 0.130746 -26.1193 0 0 0 1 

		   char *pcName = pch;
		   char *pcTrans = strstr( pch, "=" );
		   pcTrans[-1] = 0; // End file name string pcName
		   //pcTrans++; // Increment to just after equal sign

		   string filename = dirName + outputCode + pcName + ".png";// + pcTrans;
		   trainFilenames.push_back( filename );

		   char *pch = pcTrans;

			for( int j =0; j < 12; j++ )
			{
				pch = strchr( pch + 1, ' ' );
				if( !pch )
					return;
				vfTrans[j] = atof( pch );
				pch++;
			}
			trainImageData.push_back( vfTrans );
	   }
	   if( strstr( pch, "ElementDataFile = LOCAL" ) )
	   {
		   // Done reading
		   break;
	   }
	}

	int iPosition = file.tellg();
	file.close();

	FILE *infile = fopen( filename.c_str(), "rb" );
	fseek( infile, iPosition, SEEK_SET );

	unsigned char *pucImgData = new unsigned char[iImgRows*iImgCols];
	Mat mtImg(iImgRows, iImgCols, CV_8UC1, pucImgData);

	// Read & write images
	for( int i = 0; i < iImgCount; i++ )
	{
		fread( mtImg.data, 1, iImgRows*iImgCols, infile );
		//imwrite( trainFilenames[i], trainImages[i] );
		imwrite( trainFilenames[i], mtImg );
    }

	delete [] pucImgData;

	fclose( infile );
    
}

//
// convertMNCToMHA_mha()
//
// Store raw 2D mnc files into an mha file.
//
// Usage:
//
//  convertMNCToMHA_mha(
//		"C:\group1\02\pre\sweep_2a\2d\2a.txt", // Input file listing mnc file name + transforms
//		"C:\group1\02\pre\sweep_2a\2d\2a.mha", // Output file
//		"C:\group1\02\pre\sweep_2a\2d"  // Directory name for input mnc image files
//	);
//
void
convertMNCToMHA_mha(
							 const string& filename,
							 const string& outputCode
							 )
{
#ifdef WIN32
    const char dlmtr = '\\';
#else
    const char dlmtr = '/';
#endif

	vector<string> trainFilenames;
	vector< vector<float> > trainImageData;

    trainFilenames.clear();

    ifstream file( filename.c_str() );
    if ( !file.is_open() )
        return;

    size_t pos = filename.rfind(dlmtr);
    string dirName = pos == string::npos ? "" : filename.substr(0, pos) + dlmtr;

	int iImgCols = -1;
	int iImgRows = -1;
	int iImgCount = -1;

	int iCurrImage = 0;

	// Vector for reading in transforms
	vector< float > vfTrans;
	vfTrans.resize(12);

	// Read until get dimensions
    while( !file.eof() )
    {
       string str; getline( file, str );
       if( str.empty() ) break;
	   char *pch = &(str[0]);
	   if( !pch )
		   return;

	   if( strstr( pch, "DimSize =" ) )
	   {
		   if( sscanf( pch, "DimSize = %d %d %d", &iImgCols, &iImgRows, &iImgCount ) != 3 )
		   {
			   printf( "Error: could not read dimensions\n" );
			   return;
		   }
	   }
	   if( strstr( pch, "ProbeToTrackerTransform =" )
		   || strstr( pch, "UltrasoundToTrackerTransform =" ) )
	   {
		   // Parse name and transform
		   // Seq_Frame0000_ProbeToTrackerTransform = -0.224009 -0.529064 0.818481 212.75 0.52031 0.6452 0.559459 -14.0417 -0.824074 0.551188 0.130746 -26.1193 0 0 0 1 

		   char *pcName = pch;
		   char *pcTrans = strstr( pch, "=" );
		   pcTrans[-1] = 0; // End file name string pcName
		   //pcTrans++; // Increment to just after equal sign

		   string filename = dirName + outputCode + pcName + ".png";// + pcTrans;
		   trainFilenames.push_back( filename );

		   char *pch = pcTrans;

			for( int j =0; j < 12; j++ )
			{
				pch = strchr( pch + 1, ' ' );
				if( !pch )
					return;
				vfTrans[j] = atof( pch );
				pch++;
			}
			trainImageData.push_back( vfTrans );
	   }
	   if( strstr( pch, "ElementDataFile = LOCAL" ) )
	   {
		   // Done reading
		   break;
	   }
	}

	int iPosition = file.tellg();
	file.close();

	FILE *infile = fopen( filename.c_str(), "rb" );
	fseek( infile, iPosition, SEEK_SET );

	unsigned char *pucImgData = new unsigned char[iImgRows*iImgCols];
	Mat mtImg(iImgRows, iImgCols, CV_8UC1, pucImgData);

	// Read & write images
	for( int i = 0; i < iImgCount; i++ )
	{
		fread( mtImg.data, 1, iImgRows*iImgCols, infile );
		//imwrite( trainFilenames[i], trainImages[i] );
		imwrite( trainFilenames[i], mtImg );
    }

	delete [] pucImgData;

	fclose( infile );
    
}

//
// convertPngTextToMHA_mha()
//
// Convert a list text file with png images to an mha file.
//
bool
convertPngTextToMHA_mha(
				   const string& trainFilename,
				   const string& outputFilename,
				   const string& imageDirName
				   )
{
#ifdef WIN32
    const char dlmtr = '\\';
#else
    const char dlmtr = '/';
#endif

	vector<string> trainImageNames;
	vector< vector<float> > trainImageData;

    string trainDirName;
    readTrainFilenames( trainFilename, trainDirName, trainImageNames );
    if( trainImageNames.empty() )
    {
        cout << "Train image filenames can not be read." << endl << ">" << endl;
        return false;
    }

	// Fixed for MNI BYTE data
	int iImgCols = 640;
	int iImgRows = 480;

	char pcHead1[] =
	"ObjectType = Image\n \
NDims = 3\n\
AnatomicalOrientation = RAI\n\
BinaryData = True\n\
BinaryDataByteOrderMSB = False\n\
CenterOfRotation = 0 0 0\n\
CompressedData = False\n";

	// "DimSize = 640 480 13\n";

	char pcHead2 [] = 
	"ElementSpacing = 1 1 1\n\
ElementType = MET_UCHAR\n\
Offset = 0 0 0\n\
TransformMatrix = 1 0 0 0 1 0 0 0 1\n\
UltrasoundImageOrientation = MF\n\
UltrasoundImageType = BRIGHTNESS\n";

	// Output minimal MHA header information

	FILE *outfile = fopen( outputFilename.c_str(), "wb" );
	if( !outfile )
	{
		return false;
	}
	fprintf( outfile, "%s", pcHead1 );
	fprintf( outfile, "DimSize = %d %d %d\n", iImgCols, iImgRows, trainImageNames.size() );
	fprintf( outfile, "%s", pcHead2 );

	// Go through and output f
    for( size_t i = 0; i < trainImageNames.size(); i++ )
    {
		// Look for first tab after file name
		char *pch = &(trainImageNames[i][0]);
		if( !pch ) return false;

		// Read past name
		pch = strchr( pch, '\t' );
		if( !pch ) return false;
		pch++;

		// Read past index
		pch = strchr( pch, '\t' );
		if( !pch ) return false;
		pch++;

		fprintf( outfile, "Seq_Frame%4.4d_UltrasoundToTrackerTransform = %s\n", i, pch );
		fprintf( outfile, "Seq_Frame%4.4d_UltrasoundToTrackerTransformStatus = OK\n", i );
		fprintf( outfile, "Seq_Frame%4.4d_Timestamp = %d\n", i , i );
		fprintf( outfile, "Seq_Frame%4.4d_ImageStatus = OK\n", i );
    }

	fprintf( outfile, "ElementDataFile = LOCAL\n" );

	unsigned char *pucImgData = new unsigned char[iImgCols*iImgRows];

	for( size_t i = 0; i < trainImageNames.size(); i++ )
    {
		// Look for first tab after file name
		char *pchName = &(trainImageNames[i][0]);
		if( !pchName ) return false;

		// Read past name
		char *pchTab = strchr( pchName, '\t' );
		if( !pchTab ) return false;
		*pchTab = 0;
		
		// Find slash between directory and file names
		char *pchSlash = strrchr( pchName, dlmtr );
		if( !pchSlash )
			pchSlash = pchName;

		// Find extension, convert to .mnc (raw image data)
		char *pchExt = strrchr( pchSlash, '.' );
		if( !pchExt )
		{
			cout << "Error: no image name extension." << endl;
			return false;
		}
		pchExt[1] = 'm';
		pchExt[2] = 'n';
		pchExt[3] = 'c';

		// Make new name
		string filename = imageDirName + pchSlash;
		// Read in
		FILE *fileImgIn = fopen( filename.c_str(), "rb" );
		if( !fileImgIn )
		{
            cout << "Train image " << fileImgIn << " can not be read." << endl;
			return false;
		}

		// Read in raw unsigned char data from end of file
		fseek( fileImgIn, -iImgCols*iImgRows, SEEK_END );
		fread( pucImgData, 1, iImgCols*iImgRows, fileImgIn );
		fclose( fileImgIn );

		// Print out
		fwrite( pucImgData, 1, iImgCols*iImgRows, outfile );

		// Test: this works
		//FILE *outjunk = fopen( "junk.bin", "wb" );
		//fwrite( pucImgData, 1, iImgCols*iImgRows, outjunk );
		//fclose( outjunk );
	}

	delete [] pucImgData;

	fclose( outfile );

	return true;
}


bool createDetectorDescriptorMatcher( const string& detectorType, const string& descriptorType, const string& matcherType,
                                      Ptr<FeatureDetector>& featureDetector,
                                      Ptr<DescriptorExtractor>& descriptorExtractor,
                                      Ptr<DescriptorMatcher>& descriptorMatcher )
{
    cout << "< Creating feature detector..." << endl;
	//cout << " " << detectorType <<" " << featureDetector << endl; 
    featureDetector = FeatureDetector::create( detectorType );
    cout << "descriptor extractor..." << endl;
    descriptorExtractor = DescriptorExtractor::create( descriptorType );
    cout << "descriptor matcher ..." << endl;
    descriptorMatcher = DescriptorMatcher::create( matcherType );
    cout << "done >" << endl;

    bool isCreated = !( featureDetector.empty() || descriptorExtractor.empty() || descriptorMatcher.empty() );
    if( !isCreated )
        cout << "Can not create feature detector or descriptor extractor or descriptor matcher of given types." << endl << ">" << endl;

    return isCreated;
}

bool readImages( const string& queryImageName, const string& trainFilename,
                 Mat& queryImage, vector <Mat>& trainImages, vector<string>& trainImageNames )
{
    cout << "< Reading the images..." << endl;
    queryImage = imread( queryImageName, CV_LOAD_IMAGE_GRAYSCALE);
    if( queryImage.empty() )
    {
        cout << "Query image can not be read." << endl << ">" << endl;
        return false;
    }
    string trainDirName;
    readTrainFilenames( trainFilename, trainDirName, trainImageNames );
    if( trainImageNames.empty() )
    {
        cout << "Train image filenames can not be read." << endl << ">" << endl;
        return false;
    }
    int readImageCount = 0;
    for( size_t i = 0; i < trainImageNames.size(); i++ )
    {
        string filename = trainDirName + trainImageNames[i];
        Mat img = imread( filename, CV_LOAD_IMAGE_GRAYSCALE );
        if( img.empty() )
            cout << "Train image " << filename << " can not be read." << endl;
        else
            readImageCount++;
        trainImages.push_back( img );
    }
    if( !readImageCount )
    {
        cout << "All train images can not be read." << endl << ">" << endl;
        return false;
    }
    else
        cout << readImageCount << " train images were read." << endl;
    cout << ">" << endl;

    return true;
}

bool
readImageList(
				   const string& trainFilename,
				   vector <Mat>& trainImages,
				   vector<string>& trainImageNames,
				   vector< vector<float> >& trainImageData
				   )
{
    string trainDirName;
    readTrainFilenames( trainFilename, trainDirName, trainImageNames );
    if( trainImageNames.empty() )
    {
        cout << "Train image filenames can not be read." << endl << ">" << endl;
        return false;
    }
    int readImageCount = 0;
	trainImageData.resize( trainImageNames.size() );
    for( size_t i = 0; i < trainImageNames.size(); i++ )
    {
		// Read in affine transform parameters
		trainImageData[i].resize( 12 );
		// Look for first tab after file name
		char *pch = &(trainImageNames[i][0]);
		if( !pch ) return false;

		// Read past name
		pch = strchr( pch, '\t' );
		if( !pch ) return false;
		pch++;

		// Read past index
		pch = strchr( pch, '\t' );
		if( !pch ) return false;
		pch++;
		for( int j =0 ;j < 12; j++ )
		{
			while( pch[0] == ' ' || pch[0] == '\t' )
				pch++;
			//pch = strchr( pch, '\t' );
			//if( !pch ) return false;
			//pch++;
			trainImageData[i][j] = atof( pch );		
			while( pch[0] != ' ' && pch[0] != '\t' )
				pch++;

		}
		pch = &(trainImageNames[i][0]);
		pch = strchr( pch, '\t' );
		if( pch ) *pch = 0;
		pch = &(trainImageNames[i][0]);
		pch = strstr( pch, ".crop.png" );
		if( pch ) *pch = 0;

        string filename = trainDirName + trainImageNames[i];


		Mat img = imread( filename, CV_LOAD_IMAGE_GRAYSCALE );
        if( img.empty() )
            cout << "Train image " << filename << " can not be read." << endl;
        else
            readImageCount++;
        trainImages.push_back( img );
    }
    if( !readImageCount )
    {
        cout << "All train images can not be read." << endl << ">" << endl;
        return false;
    }
    else
        cout << readImageCount << " train images were read." << endl;
    cout << ">" << endl;

    return true;
}

bool
readImageList_mha(
				   const string& trainFilename,
				   vector <Mat>& trainImages,
				   vector<string>& trainImageNames,
				   vector< vector<float> >& trainImageData
				   )
{
    string trainDirName;
    readTrainTransforms_mha( trainFilename, trainDirName, trainImageNames, trainImages, trainImageData );
	int readImageCount = trainImages.size();

	if( !readImageCount )
    {
        cout << "All train images can not be read." << endl << ">" << endl;
        return false;
    }
    else
        cout << readImageCount << " train images were read." << endl;
    cout << ">" << endl;

    return true;
}



bool
readImageListDetectCompute(
								const string& trainFilename,
								vector <Mat>& trainImages,
								vector<string>& trainImageNames,
								Ptr<FeatureDetector>& featureDetector,
								Ptr<DescriptorExtractor>& descriptorExtractor,
								vector<vector<KeyPoint> > &trainKeypoints,
								vector<Mat> &trainDescriptors,
								Mat &maskImage
								)
{
    string trainDirName;
    readTrainFilenames( trainFilename, trainDirName, trainImageNames );
    if( trainImageNames.empty() )
    {
        cout << "Train image filenames can not be read." << endl << ">" << endl;
        return false;
    }
    int readImageCount = 0;

	trainKeypoints.resize( trainImageNames.size() );
    for( size_t i = 0; i < trainImageNames.size(); i++ )
    {
        string filename = trainDirName + trainImageNames[i];
        Mat img = imread( filename, CV_LOAD_IMAGE_GRAYSCALE );
        if( img.empty() )
            cout << "Train image " << filename << " can not be read." << endl;
        else
            readImageCount++;
        trainImages.push_back( img );

		char *pch = strchr( &filename[0], '.' );
		if( !pch )
		{
			printf( "Error: strange file name, no image extension\n" );
			return false;
		}

		featureDetector->detect( img, trainKeypoints[i], maskImage );
    }
	descriptorExtractor->compute( trainImages, trainKeypoints, trainDescriptors );

    if( !readImageCount )
    {
        cout << "All train images can not be read." << endl << ">" << endl;
        return false;
    }
    else
        cout << readImageCount << " train images were read." << endl;
    cout << ">" << endl;

    return true;
}

void detectKeypoints( const Mat& queryImage, vector<KeyPoint>& queryKeypoints,
                      const vector<Mat>& trainImages, vector<vector<KeyPoint> >& trainKeypoints,
                      Ptr<FeatureDetector>& featureDetector )
{
    cout << endl << "< Extracting keypoints from images..." << endl;
    featureDetector->detect( queryImage, queryKeypoints );
    featureDetector->detect( trainImages, trainKeypoints );
    cout << ">" << endl;
}

void computeDescriptors( const Mat& queryImage, vector<KeyPoint>& queryKeypoints, Mat& queryDescriptors,
                         const vector<Mat>& trainImages, vector<vector<KeyPoint> >& trainKeypoints, vector<Mat>& trainDescriptors,
                         Ptr<DescriptorExtractor>& descriptorExtractor )
{
    cout << "< Computing descriptors for keypoints..." << endl;
    descriptorExtractor->compute( queryImage, queryKeypoints, queryDescriptors );
    descriptorExtractor->compute( trainImages, trainKeypoints, trainDescriptors );
    cout << ">" << endl;
}

void matchDescriptors( const Mat& queryDescriptors, vector<DMatch>& matches, Ptr<DescriptorMatcher>& descriptorMatcher )
{
	// Assumes training descriptors have already been added to descriptorMatcher
    cout << "< Set train descriptors collection in the matcher and match query descriptors to them..." << endl;
    //descriptorMatcher->add( trainDescriptors );
    descriptorMatcher->match( queryDescriptors, matches );
    CV_Assert( queryDescriptors.rows == (int)matches.size() || matches.empty() );
    cout << ">" << endl;
}

void matchDescriptorsKNN( const Mat& queryDescriptors, vector<vector<DMatch> >& matches, Ptr<DescriptorMatcher>& descriptorMatcher, int k)
{
	// Assumes training descriptors have already been added to descriptorMatcher
    cout << "< Set train descriptors collection in the matcher and match query descriptors to them..." << endl;
    //descriptorMatcher->add( trainDescriptors );
	descriptorMatcher->knnMatch( queryDescriptors, matches, k );
    CV_Assert( queryDescriptors.rows == (int)matches.size() || matches.empty() );
    cout << ">" << endl;
}

void saveResultImage( const Mat& queryImage, const vector<KeyPoint>& queryKeypoints,
                       const Mat& trainImage, const vector<KeyPoint>& trainKeypoints,
                       const vector<DMatch>& matches,
					   int iTrainImageIndex,
					   const string& outputImageName,
					   const string& resultDir
					   )
{
    cout << "< Save results..." << endl;
    Mat drawImg;
    vector<char> mask;
    mask.resize( matches.size() );
    fill( mask.begin(), mask.end(), 0 );
    for( size_t i = 0; i < matches.size(); i++ )
    {
        if( abs( matches[i].imgIdx - iTrainImageIndex ) <= 1 )
            mask[i] = 1;
    }

	maskMatchesByTrainImgIdx( matches, iTrainImageIndex, mask );
	maskMatchesByMatchPresent( matches, mask );
	drawMatches( queryImage, queryKeypoints, trainImage, trainKeypoints,
		matches, drawImg, Scalar::all(-1), Scalar::all(-1), mask, DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS | DrawMatchesFlags::DRAW_RICH_KEYPOINTS );
	int iLast = outputImageName.find_last_of( "\\" );
	int iSize  = outputImageName.size();
	string fname = outputImageName.substr( iLast + 1, iSize-iLast-1);
	string filename = resultDir + "\\single_res_" + fname;
	imwrite( filename, drawImg );
}

void
saveResultData(
			FILE *outfile1,
			FILE *outfile2,
			vector<string>& queryImagesNames,
			vector<vector<KeyPoint> >& queryKeypoints,
			vector<string>& trainImagesNames,
			vector<vector<KeyPoint> >& trainKeypoints,
			vector<DMatch>& matches,
			int iQueryImageIndex
			)
{
	string& queryName = queryImagesNames[iQueryImageIndex];
	vector<KeyPoint>& queryKeys = queryKeypoints[iQueryImageIndex];

    for( int i = 0; i < matches.size(); i++ )
    {
		if( matches[i].imgIdx >= 0 )
		{
			KeyPoint &queryKey = queryKeys[matches[i].queryIdx];
			KeyPoint &trainKey = trainKeypoints[matches[i].imgIdx][matches[i].trainIdx];

			string& trainName = trainImagesNames[matches[i].imgIdx];

			fprintf( outfile1, "%d|%s 0 \%f \%f|%s 0 %f %f\n",
				matches.size(),
				queryName.c_str(), queryKey.pt.y, queryKey.pt.x,
				trainName.c_str(), trainKey.pt.y, trainKey.pt.x );
		}
	}
	return;
}

void saveResultImages( const Mat& queryImage, const vector<KeyPoint>& queryKeypoints,
                       const vector<Mat>& trainImages, const vector<vector<KeyPoint> >& trainKeypoints,
                       const vector<DMatch>& matches, const vector<string>& trainImagesNames, const string& resultDir,
					   int iIndex = 0 )
{
    cout << "< Save results..." << endl;
    Mat drawImg;
    vector<char> mask;
	int iMatches;
    for( size_t i = 0; i < trainImages.size(); i++ )
    {
        if( !trainImages[i].empty() )
        {
            iMatches = maskMatchesByTrainImgIdx( matches, i, mask );
			iMatches = maskMatchesByMatchPresent( matches, mask );
			if( iMatches <= 0 )
				continue;
            drawMatches( queryImage, queryKeypoints, trainImages[i], trainKeypoints[i],
				matches, drawImg, Scalar::all(-1), Scalar::all(-1), mask, DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS | DrawMatchesFlags::DRAW_RICH_KEYPOINTS );
			int iLast = trainImagesNames[i].find_last_of( "\\" );
			int iSize  = trainImagesNames[i].size();
			string fname = trainImagesNames[i].substr( iLast + 1, iSize-iLast-1);
			char pcIndex[10];
			sprintf( pcIndex, "%3.3d_", iIndex );
            string filename = resultDir + "\\res_" + pcIndex + fname;
            if( !imwrite( filename, drawImg ) )
                cout << "Image " << filename << " can not be saved (may be because directory " << resultDir << " does not exist)." << endl;

//			Mat outImg1, outImg2;
//			outImg1 = drawImg( Rect(0, 0, queryImage.cols, queryImage.rows) );
//			outImg2 = drawImg( Rect(queryImage.cols, 0, trainImages[i].cols, trainImages[i].rows) );

   //         string filename1 = resultDir + "\\res_1_" + fname + "_1.png";
			//imwrite( filename1, outImg1 );
   //         string filename2 = resultDir + "\\res_2_" + fname + "_2.png";
			//imwrite( filename2, outImg2 );


        }
    }
    cout << ">" << endl;
}

//
// saveResultImages()
//
// Output multiple matches.
//
void saveResultImages( const Mat& queryImage, const vector<KeyPoint>& queryKeypoints,
                       const vector<Mat>& trainImages, const vector<vector<KeyPoint> >& trainKeypoints,
                       const vector< vector<DMatch> >& mmatches, const vector<string>& trainImagesNames, const string& resultDir )
{
    cout << "< Save results..." << endl;

	vector<DMatch> matches;
	matches.resize( mmatches.size() );
	for( int i = 0; i < mmatches.size(); i++ )
	{
		matches[i] = mmatches[i][0];
	}

		saveResultImages( queryImage, queryKeypoints, trainImages, trainKeypoints,
			matches, trainImagesNames, resultDir );

   // Mat drawImg;
   // vector<char> mask;
   // for( size_t i = 0; i < trainImages.size(); i++ )
   // {
   //     if( !trainImages[i].empty() )
   //     {
   //         //maskMatchesByTrainImgIdx( matches, i, mask );
   //         drawMatches( queryImage, queryKeypoints, trainImages[i], trainKeypoints[i],
   //                      matches, drawImg, Scalar::all(-1), Scalar::all(-1), mask );
			//int iLast = trainImagesNames[i].find_last_of( "\\" );
			//int iSize  = trainImagesNames[i].size();
			//string fname = trainImagesNames[i].substr( iLast + 1, iSize-iLast-1);
   //         string filename = resultDir + "\\res_" + fname;
   //         if( !imwrite( filename, drawImg ) )
   //             cout << "Image " << filename << " can not be saved (may be because directory " << resultDir << " does not exist)." << endl;
   //     }
   // }
    cout << ">" << endl;
}

//
// applyImageMask()
//
// Apply binary mask.
//
void
applyImageMask(
			   const Mat & maskImage,
			   vector<Mat>& vecImages
		  )
{
	for( int i = 0; i < vecImages.size(); i++ )
	{
		for( int y = 0; y < vecImages[i].rows; y++ )
		{
			for( int x = 0; x < vecImages[i].cols; x++ )
			{
				if( maskImage.data[ y*maskImage.step + x ] <= 0 )
				{
					vecImages[i].data[ y*vecImages[i].step + x ] = 0;
				}
			}
		}
	}
}

int
_readDescriptorDatabase(
						   char *pcFileWithImages,
						vector<Mat> &trainImages,
						vector<string> &trainImagesNames,
						vector< vector<float> > &trainImagesTransform,
						vector<vector<KeyPoint> > &trainKeypoints,
						vector<Mat> &trainDescriptors
						  )
{
    string fileWithTrainImages = pcFileWithImages;

	char *pch = &(fileWithTrainImages[0]);
	if( !pch ) return false;
	if( strstr( pch, ".mha" ) )
	{
		if( !readImageList_mha( fileWithTrainImages, trainImages, trainImagesNames, trainImagesTransform ) )
		{
			printf( "Error: could not read training images\n" );
			return -1;
		}
	}
	else if( !readImageList( fileWithTrainImages, trainImages, trainImagesNames, trainImagesTransform ) )
	{
		printf( "Error: could not read training images\n" );
		return -1;
	}

	pch = strstr( &fileWithTrainImages[0], ".txt" );
	if( !pch )
	{
		pch = strstr( &fileWithTrainImages[0], ".mha" );

		if( !pch )
		{
			printf( "Error: strange file name, no extension .txt\n" );
			return false;
		}
	}

	fileWithTrainImages += "extra extra extra space";
	pch = strrchr( &fileWithTrainImages[0], '.' );

	FILE *outfile;
	int iValue;


	sprintf( pch, "%s", ".surf.key" );
	outfile = fopen( &fileWithTrainImages[0], "rb" );
	if( !outfile )
	{
		printf( "Error: no key file found\n" );
		return false;
	}
	fread( &iValue, sizeof(iValue), 1, outfile );	
	assert( iValue == trainImages.size() );
	trainKeypoints.resize( iValue ); // Number of vector<Keypoint>/Images
	for( int i = 0; i < trainKeypoints.size(); i++ )
	{
		fread( &iValue, sizeof(iValue), 1, outfile );
		trainKeypoints[i].resize(iValue); // Size of vector<Keypoint>
		for( int j = 0; j < trainKeypoints[i].size(); j++ )
		{
			fread( &(trainKeypoints[i][j]), sizeof(trainKeypoints[i][j]), 1, outfile );
		}
	}
	fclose( outfile );

	int iReturn;
	sprintf( pch, "%s", ".surf.dsc" );
	outfile = fopen( &fileWithTrainImages[0], "rb" );
	fread( &iValue, sizeof(iValue), 1, outfile );
	assert( iValue == trainImages.size() );
	trainDescriptors.resize(iValue); // Number of Mat/Images
	for( int i = 0; i < trainDescriptors.size(); i++ )
	{
		int iRows;
		int iCols;
		int iType;
		fread( &iRows, sizeof(iRows), 1, outfile );
		fread( &iCols, sizeof(iCols), 1, outfile );
		fread( &iType, sizeof(iType), 1, outfile );
		if( iRows > 0 && iCols > 0 )
		{
			assert( iType == CV_32FC1 );
			trainDescriptors[i].create( iRows, iCols, iType );
			iReturn = fread( trainDescriptors[i].data, sizeof(float), iRows*iCols, outfile );
		}
		else
		{
			printf( "Zero feature image: %d\n", i );
		}
	}
	fclose( outfile );

	return 0;
}

//
// _readDescriptorDatabase_range()
//
// Read in features/descriptors within a start/end range.
// The idea is to reduce input descriptors for a more effective correspondence / search
//
int
_readDescriptorDatabase_range(
						   char *pcFileWithImages,
						vector<Mat> &trainImages,
						vector<string> &trainImagesNames,
						vector< vector<float> > &trainImagesTransform,
						vector<vector<KeyPoint> > &trainKeypoints,
						vector<Mat> &trainDescriptors,

						// Option to add start / end range 
						int iStart=0,
						int iEnd=-1
						  )
{
    string fileWithTrainImages = pcFileWithImages;

	char *pch = &(fileWithTrainImages[0]);
	if( !pch ) return false;
	if( strstr( pch, ".mha" ) )
	{
		if( !readImageList_mha( fileWithTrainImages, trainImages, trainImagesNames, trainImagesTransform ) )
		{
			printf( "Error: could not read training images\n" );
			return -1;
		}
	}
	else if( !readImageList( fileWithTrainImages, trainImages, trainImagesNames, trainImagesTransform ) )
	{
		printf( "Error: could not read training images\n" );
		return -1;
	}

	pch = strstr( &fileWithTrainImages[0], ".txt" );
	if( !pch )
	{
		pch = strstr( &fileWithTrainImages[0], ".mha" );

		if( !pch )
		{
			printf( "Error: strange file name, no extension .txt\n" );
			return false;
		}
	}

	fileWithTrainImages += "extra extra extra space";
	pch = strrchr( &fileWithTrainImages[0], '.' );

	FILE *outfile;
	int iValue;


	sprintf( pch, "%s", ".surf.key" );
	outfile = fopen( &fileWithTrainImages[0], "rb" );
	if( !outfile )
	{
		printf( "Error: no key file found\n" );
		return false;
	}
	fread( &iValue, sizeof(iValue), 1, outfile );	
	assert( iValue == trainImages.size() );
	trainKeypoints.resize( iValue ); // Number of vector<Keypoint>/Images
	for( int i = 0; i < trainKeypoints.size(); i++ )
	{
		fread( &iValue, sizeof(iValue), 1, outfile );
		trainKeypoints[i].resize(iValue); // Size of vector<Keypoint>
		for( int j = 0; j < trainKeypoints[i].size(); j++ )
		{
			fread( &(trainKeypoints[i][j]), sizeof(trainKeypoints[i][j]), 1, outfile );
		}
	}
	fclose( outfile );

	int iReturn;
	sprintf( pch, "%s", ".surf.dsc" );
	outfile = fopen( &fileWithTrainImages[0], "rb" );
	fread( &iValue, sizeof(iValue), 1, outfile );
	assert( iValue == trainImages.size() );
	trainDescriptors.resize(iValue); // Number of Mat/Images
	for( int i = 0; i < trainDescriptors.size(); i++ )
	{
		int iRows;
		int iCols;
		int iType;
		fread( &iRows, sizeof(iRows), 1, outfile );
		fread( &iCols, sizeof(iCols), 1, outfile );
		fread( &iType, sizeof(iType), 1, outfile );
		if( iRows > 0 && iCols > 0 )
		{
			assert( iType == CV_32FC1 );
			trainDescriptors[i].create( iRows, iCols, iType );
			iReturn = fread( trainDescriptors[i].data, sizeof(float), iRows*iCols, outfile );
		}
		else
		{
			printf( "Zero feature image: %d\n", i );
		}
	}
	fclose( outfile );

	// Here, erase selectively

	if( iEnd == -1 || iEnd <= iStart )
	{
		iEnd = trainImages.size();
	}
	
	int iWidth = iEnd-iStart;
	
	if( iWidth > 0 )
	{
		trainImages.erase( trainImages.begin()+iEnd, trainImages.end() );
		trainImagesNames.erase( trainImagesNames.begin()+iEnd, trainImagesNames.end() );
		trainImagesTransform.erase( trainImagesTransform.begin()+iEnd, trainImagesTransform.end() );
		trainKeypoints.erase( trainKeypoints.begin()+iEnd, trainKeypoints.end() );
		trainDescriptors.erase( trainDescriptors.begin()+iEnd, trainDescriptors.end() );
	}

	if( iStart > 0 && iStart < trainImages.size() )
	{
		trainImages.erase( trainImages.begin(), trainImages.begin()+iStart );
		trainImagesNames.erase( trainImagesNames.begin(), trainImagesNames.begin()+iStart );
		trainImagesTransform.erase( trainImagesTransform.begin(), trainImagesTransform.begin()+iStart );
		trainKeypoints.erase( trainKeypoints.begin(), trainKeypoints.begin()+iStart );
		trainDescriptors.erase( trainDescriptors.begin(), trainDescriptors.begin()+iStart );
	}


	return 0;
}



//
// main_prepareDescriptorDatabase()
//
// Produce & save descriptor database.
//
int
main_prepareDescriptorDatabase(
							   char *pcDetectorType,
							   char *pcDescriptorType,
							   char *pcFileWithImages,
							   char *pcMaskImage
						  )
{
	string detectorType = defaultDetectorType;
    string descriptorType = defaultDescriptorType;
    string matcherType = defaultMatcherType;
    string fileWithQueryImages = defaultQueryImageName;
    string fileWithTrainImages = defaultFileWithTrainImages;
    string fileWithBogusImages = "";
    string fileMaskImage = "";
    string dirToSaveResImages = defaultDirToSaveResImages;

	detectorType = pcDetectorType;
    descriptorType = pcDescriptorType;
    fileWithTrainImages = pcFileWithImages;
	fileMaskImage = pcMaskImage;

    Ptr<FeatureDetector> featureDetector;
    Ptr<DescriptorExtractor> descriptorExtractor;
    Ptr<DescriptorMatcher> descriptorMatcher;
    Ptr<DescriptorMatcher> descriptorMatcherBogus;
 
	printf( "Creating detector/descriptor/matcher: %s %s %s\n", detectorType.data(), descriptorType.data(), matcherType.data() );
	if( !createDetectorDescriptorMatcher( detectorType, descriptorType, matcherType, featureDetector, descriptorExtractor, descriptorMatcher ) )
    {
		printf( "Error: could not create descriptor extractor matcher\n" );
        return false;
    }

	printf( "Reading mask: %s\n", fileMaskImage.data() );
    Mat maskImage = imread( fileMaskImage, CV_LOAD_IMAGE_GRAYSCALE);
	if( maskImage.empty() )
	{
		cout << "Mask image can not be read:" << fileMaskImage << endl << ">" << endl;
        return false;
	}

    vector<Mat> trainImages;
    vector<string> trainImagesNames;
    vector<vector<float>> trainTransforms;
    vector<vector<KeyPoint> > trainKeypoints;
    vector<Mat> trainDescriptors;

	printf( "Reading image list: %s\n", fileWithTrainImages.data() );

	// Look for first tab after file name
	char *pch = &(fileWithTrainImages[0]);
	if( !pch ) return false;
	if( strstr( pch, ".mha" ) )
	{
		if( !readImageList_mha( fileWithTrainImages, trainImages, trainImagesNames, trainTransforms ) )
		{
			printf( "Error: could not read training images\n" );
			return -1;
		}
	}
	else if( !readImageList( fileWithTrainImages, trainImages, trainImagesNames, trainTransforms ) )
	{
		printf( "Error: could not read training images\n" );
        return -1;
    }

	trainKeypoints.resize( trainImages.size() );
	for( int i = 0; i < trainImages.size(); i++ )
	{
		printf( "%d\n", i );
		featureDetector->detect( trainImages[i], trainKeypoints[i], maskImage );
	}
	descriptorExtractor->compute( trainImages, trainKeypoints, trainDescriptors );

	// Original text file input (AMIGO)
	pch = strstr( &fileWithTrainImages[0], ".txt" );
	if( !pch )
	{
		// PLUS mha file input
		pch = strstr( &fileWithTrainImages[0], ".mha" );
	}
	if( !pch )
	{
		printf( "Error: strange file name, no extension .txt\n" );
		return false;
	}
	pch = strrchr( &fileWithTrainImages[0], '.' );

	FILE *outfile;
	int iValue;

	sprintf( pch, "%s", ".surf.key" );
	outfile = fopen( &fileWithTrainImages[0], "wb" );
	iValue = trainKeypoints.size(); // Number of vector<Keypoint>/Images
	fwrite( &iValue, sizeof(iValue), 1, outfile );
	for( int i = 0; i < trainKeypoints.size(); i++ )
	{
		int iValue = trainKeypoints[i].size(); // Size of vector<Keypoint>
		fwrite( &iValue, sizeof(iValue), 1, outfile );
		for( int j = 0; j < trainKeypoints[i].size(); j++ )
		{
			fwrite( &(trainKeypoints[i][j]), sizeof(trainKeypoints[i][j]), 1, outfile );
		}
	}
	fclose( outfile );
	
	int iReturn;
	sprintf( pch, "%s", ".surf.dsc" );
	outfile = fopen( &fileWithTrainImages[0], "wb" );
	iValue = trainDescriptors.size(); // Number of Mat/Images
	fwrite( &iValue, sizeof(iValue), 1, outfile );
	for( int i = 0; i < trainDescriptors.size(); i++ )
	{
		int iRows = trainDescriptors[i].rows; // Size of vector<Keypoint>
		int iCols = trainDescriptors[i].cols; // Size of vector<Keypoint>
		int iType = trainDescriptors[i].type(); // Size of vector<Keypoint>
		assert( iType == CV_32FC1 );
		fwrite( &iRows, sizeof(iRows), 1, outfile );
		fwrite( &iCols, sizeof(iCols), 1, outfile );
		fwrite( &iType, sizeof(iType), 1, outfile );
		iReturn = fwrite( trainDescriptors[i].data, sizeof(float), iRows*iCols, outfile );
	}
	fclose( outfile );

	
    vector<Mat> trainImages2;
    vector<string> trainImagesNames2;
    vector<vector<float>> trainTransforms2;
    vector<vector<KeyPoint> > trainKeypoints2;
    vector<Mat> trainDescriptors2;
	//main_readDescriptorDatabase( pcFileWithImages, trainImages2, trainImagesNames2, trainTransforms2, trainKeypoints2,trainDescriptors2 );

	return 0;
}


//
// addToOutput()
//
// Create linear coefficients files, to be solved via SVD, etc.
//
// Outputs are vecX and vecY
//
int
addToOutput(
			const vector<DMatch> &matches,

			vector<float>  &queryImageTransform,
			vector<KeyPoint> &queryKeypoints,

			vector< vector<float> > &trainImagesTransforms,
			vector<vector<KeyPoint> > &trainKeypoints,

			vector< vector<float> > &vecX,
			vector< vector<float> > &vecY
			)
{
	float pfMatQuery[16];
	float pfMatQueryInv[16];

	float pfMatResult[16];

	float pfMatTrain[16];
	float pfMatTrainInv[16];

	cv::Mat matQuery( 4, 4, CV_32FC1, pfMatQuery );
	cv::Mat matQueryInv( 4, 4, CV_32FC1, pfMatQueryInv );

	cv::Mat matTrain( 4, 4, CV_32FC1, pfMatTrain );
	cv::Mat matTrainInv( 4, 4, CV_32FC1, pfMatTrainInv );

	memset( pfMatQuery, 0, sizeof(pfMatQuery) );	pfMatQuery[15] = 1;
	memset( pfMatTrain, 0, sizeof(pfMatTrain) );	pfMatTrain[15] = 1;

	// Copy and invert query image transform 
	memcpy( pfMatQuery, (&(queryImageTransform[0])), sizeof(float)*12 );
	cv::invert( matQuery, matQueryInv );

	//cvMatMul( (CvArr*)(&(CvMat)matQuery), (CvArr*)(&(CvMat)matQueryInv), (CvArr*)(&(CvMat)matResult) );

	vector< float > vec_x;
	vec_x.resize(9);
	vector< float > vec_y;
	vec_y.resize(1);

	for( int i = 0; i < matches.size(); i++ )
	{
			int iIdx = matches[i].trainIdx;
			int iImg = matches[i].imgIdx;
			if( iImg >= 0 && iIdx >= 0 )
			{
				KeyPoint kp1 = queryKeypoints[matches[i].imgIdx];
				KeyPoint kp2 = trainKeypoints[iImg][iIdx];

				memcpy( pfMatTrain, (&(trainImagesTransforms[iImg][0])), sizeof(float)*12 );
				cv::invert( matTrain, matTrainInv );

				for( int row = 0; row < 3; row++ )
				{
					vec_x[0] = pfMatQueryInv[4*row+0]*kp1.pt.x - pfMatTrainInv[4*row+0]*kp2.pt.x;
					vec_x[1] = pfMatQueryInv[4*row+0]*kp1.pt.y - pfMatTrainInv[4*row+0]*kp2.pt.y;
					
					vec_x[2] = pfMatQueryInv[4*row+1]*kp1.pt.x - pfMatTrainInv[4*row+1]*kp2.pt.x;
					vec_x[3] = pfMatQueryInv[4*row+1]*kp1.pt.y - pfMatTrainInv[4*row+1]*kp2.pt.y;

					vec_x[4] = pfMatQueryInv[4*row+2]*kp1.pt.x - pfMatTrainInv[4*row+2]*kp2.pt.x;
					vec_x[5] = pfMatQueryInv[4*row+2]*kp1.pt.y - pfMatTrainInv[4*row+2]*kp2.pt.y;
					
					vec_x[6] = pfMatQueryInv[4*row+0] - pfMatTrainInv[4*row+0];
					vec_x[7] = pfMatQueryInv[4*row+1] - pfMatTrainInv[4*row+1];
					vec_x[8] = pfMatQueryInv[4*row+2] - pfMatTrainInv[4*row+2];

					vec_y[0] = pfMatTrainInv[4*row+3] - pfMatQueryInv[4*row+3];

					vecX.push_back( vec_x );
					vecY.push_back( vec_y );
				}

				// Save
			}
	}

	//// Save resulting transform
	//outfileTransform = fopen( "match_transform.txt", "a+" );

	//fprintf( outfileTransform, "A:\t" );
	//for( int k = 0; k < queryImagesTransform[i].size(); k++ )
	//{
	//	fprintf( outfileTransform, "%f\t", queryImagesTransform[i][k] );
	//}
	//fprintf( outfileTransform, "B:\t" );
	//for( int k = 0; k < trainImagesTransform[iMaxIndex].size(); k++ )
	//{
	//	fprintf( outfileTransform, "%f\t", trainImagesTransform[iMaxIndex][k] );
	//}
	//fprintf( outfileTransform, "\n" );
	//fclose( outfileTransform );

	return 0;
}

//
//
//
// addToOutputNonLinear()
//
// Output matrices and corresponding points to non-linear solver.
//
//
int
addToOutputNonLinear(
			const vector<DMatch> &matches,

			vector<float>  &queryImageTransform,
			vector<KeyPoint> &queryKeypoints,

			vector< vector<float> > &trainImagesTransforms,
			vector<vector<KeyPoint> > &trainKeypoints,

			vector< vector<float> > &vecX,
			vector< vector<float> > &vecY
			)
{
	float pfMatQuery[16];
	float pfMatQueryInv[16];

	float pfMatVec[4];
	float pfMatResult[16];

	float pfMatTrain[16];
	float pfMatTrainInv[16];

	cv::Mat matQuery( 4, 4, CV_32FC1, pfMatQuery );
	cv::Mat matTrain( 4, 4, CV_32FC1, pfMatTrain );

	cv::Mat matVec( 4, 1, CV_32FC1, pfMatVec );
	cv::Mat matRes1( 4, 1, CV_32FC1, pfMatResult+0 );
	cv::Mat matRes2( 4, 1, CV_32FC1, pfMatResult+4 );

	memset( pfMatQuery, 0, sizeof(pfMatQuery) );	pfMatQuery[15] = 1;
	memset( pfMatTrain, 0, sizeof(pfMatTrain) );	pfMatTrain[15] = 1;

	// Copy and invert query image transform 
	memcpy( pfMatQuery, (&(queryImageTransform[0])), sizeof(float)*12 );


	vector< float > vecXInner;// Store first 3 rows of matrix (12), then input point (3)
	vector< float > vecYInner;
	vecXInner.resize( 15 );
	vecYInner.resize( 15 );

	for( int i = 0; i < matches.size(); i++ )
	{
		int iIdxQuery = matches[i].queryIdx;
		int iIdxTrain = matches[i].trainIdx;
		int iImg = matches[i].imgIdx;
		if( iImg >= 0 && iIdxTrain >= 0 )
		{
			KeyPoint kp1 = queryKeypoints[iIdxQuery];
			KeyPoint kp2 = trainKeypoints[iImg][iIdxTrain];

			memcpy( pfMatTrain, (&(trainImagesTransforms[iImg][0])), sizeof(float)*12 );

			for( int row = 0; row < 3; row++ )
			{
				for( int col = 0; col < 4; col++ )
				{
					vecXInner[4*row+col] = pfMatQuery[4*row+col];
					vecYInner[4*row+col] = pfMatTrain[4*row+col];
				}
			}
			vecXInner[4*3+0] = kp1.pt.x;
			vecXInner[4*3+1] = kp1.pt.y;
			vecXInner[4*3+2] = 0;

			vecYInner[4*3+0] = kp2.pt.x;
			vecYInner[4*3+1] = kp2.pt.y;
			vecYInner[4*3+2] = 0;

			vecX.push_back( vecXInner );
			vecY.push_back( vecYInner );

			// For MR test data

			pfMatVec[0] = kp1.pt.x;
			pfMatVec[1] = kp1.pt.y;
			pfMatVec[2] = 0;
			pfMatVec[3] = 1;
			cvMatMul( (CvArr*)(&(CvMat)matQuery), (CvArr*)(&(CvMat)matVec), (CvArr*)(&(CvMat)matRes1) );
			
			pfMatVec[0] = kp2.pt.x;
			pfMatVec[1] = kp2.pt.y;
			pfMatVec[2] = 0;
			pfMatVec[3] = 1;
			cvMatMul( (CvArr*)(&(CvMat)matTrain), (CvArr*)(&(CvMat)matVec), (CvArr*)(&(CvMat)matRes2) );

		}
	}

	//// Save resulting transform
	//outfileTransform = fopen( "match_transform.txt", "a+" );

	//fprintf( outfileTransform, "A:\t" );
	//for( int k = 0; k < queryImagesTransform[i].size(); k++ )
	//{
	//	fprintf( outfileTransform, "%f\t", queryImagesTransform[i][k] );
	//}
	//fprintf( outfileTransform, "B:\t" );
	//for( int k = 0; k < trainImagesTransform[iMaxIndex].size(); k++ )
	//{
	//	fprintf( outfileTransform, "%f\t", trainImagesTransform[iMaxIndex][k] );
	//}
	//fprintf( outfileTransform, "\n" );
	//fclose( outfileTransform );

	return 0;
}


int
test_output(
			
			const vector<DMatch> &matches,

			vector<float>  &queryImageTransform,
			vector<KeyPoint> &queryKeypoints,

			vector< vector<float> > &trainImagesTransforms,
			vector<vector<KeyPoint> > &trainKeypoints
			)
{
	// This matrix obtained via SVD
	float pfBMatrix[16] = //{-0.000383,-0.005676,0,0.000770,-0.003143,0.000443,0,0.004358,-9.016907,13.886870,0,-46.222257, 0, 0, 0, 1};
	{-0.000383, -0.005676, 0, -9.016907,
	0.000770, -0.003143, 0, 13.886870,
	0.000443, 0.004358, 0, 	-46.222257,
	0, 0, 0, 1 };


	float pfMatQuery[16];
	float pfMatQueryInv[16];

	float pfMatResult[16];

	float pfMatTrain[16];
	float pfMatTrainInv[16];

	cv::Mat matBMatrix( 4, 4, CV_32FC1, pfBMatrix );

	cv::Mat matQuery( 4, 4, CV_32FC1, pfMatQuery );
	cv::Mat matQueryInv( 4, 4, CV_32FC1, pfMatQueryInv );

	cv::Mat matTrain( 4, 4, CV_32FC1, pfMatTrain );
	cv::Mat matTrainInv( 4, 4, CV_32FC1, pfMatTrainInv );

	memset( pfMatQuery, 0, sizeof(pfMatQuery) );	pfMatQuery[15] = 1;
	memset( pfMatTrain, 0, sizeof(pfMatTrain) );	pfMatTrain[15] = 1;

	// Copy and multiply query image transform 
	memcpy( pfMatQuery, (&(queryImageTransform[0])), sizeof(float)*12 );
	cvMatMul( (CvArr*)(&(CvMat)matQuery), (CvArr*)(&(CvMat)matBMatrix), (CvArr*)(&(CvMat)matQueryInv) );

	float pfMatQueryIn[4];
	float pfMatQueryOt[4];
	float pfMatTrainIn[4];
	float pfMatTrainOt[4];
	cv::Mat matQueryIn( 4, 1, CV_32FC1, pfMatQueryIn );
	cv::Mat matQueryOt( 4, 1, CV_32FC1, pfMatQueryOt );
	cv::Mat matTrainIn( 4, 1, CV_32FC1, pfMatTrainIn );
	cv::Mat matTrainOt( 4, 1, CV_32FC1, pfMatTrainOt );

	for( int i = 0; i < matches.size(); i++ )
	{
			int iIdx = matches[i].trainIdx;
			int iImg = matches[i].imgIdx;
			if( iImg >= 0 && iIdx >= 0 )
			{
				KeyPoint kp1 = queryKeypoints[i];
				KeyPoint kp2 = trainKeypoints[iImg][iIdx];

				// Copy and multiply train image transform
				memcpy( pfMatTrain, (&(trainImagesTransforms[iImg][0])), sizeof(float)*12 );
				cvMatMul( (CvArr*)(&(CvMat)matTrain), (CvArr*)(&(CvMat)matBMatrix), (CvArr*)(&(CvMat)matTrainInv) );

				pfMatQueryIn[0] = kp1.pt.x;
				pfMatQueryIn[1] = kp1.pt.y;
				pfMatQueryIn[2] = 0;
				pfMatQueryIn[3] = 1;
				cvMatMul( (CvArr*)(&(CvMat)matQueryInv), (CvArr*)(&(CvMat)matQueryIn), (CvArr*)(&(CvMat)matQueryOt) );

				pfMatTrainIn[0] = kp2.pt.x;
				pfMatTrainIn[1] = kp2.pt.y;
				pfMatTrainIn[2] = 0;
				pfMatTrainIn[3] = 1;
				cvMatMul( (CvArr*)(&(CvMat)matTrainInv), (CvArr*)(&(CvMat)matTrainIn), (CvArr*)(&(CvMat)matTrainOt) );

			}
	}

	return 1;
}

int
writeArray(
		  char *pcFileNameBase,
		  vector< vector< float > > &vvData,
		  char *pcFormat
		  )
{
	FILE *infile;
	char pcFileName[300];
	int iReturn;

	// Open header
	sprintf( pcFileName, "%s.txt", pcFileNameBase );
	infile = fopen( pcFileName, "wt" );
	if( !infile )
	{
		return -1;
	}
	int iRows = vvData.size();
	int iCols = vvData[0].size();
	iReturn = fprintf( infile, "cols:\t%d\n", vvData[0].size() );
	iReturn = fprintf( infile, "rows:\t%d\n", vvData.size() );
	iReturn = fprintf( infile, "format:\t%s\n", pcFormat );
	fclose( infile );

	// Open data
	sprintf( pcFileName, "%s.bin", pcFileNameBase );
	infile = fopen( pcFileName, "wb" );
	if( !infile )
	{
		return -1;
	}

	if( pcFormat[0] == 'i' )
	{
		for( int i = 0; i < iRows; i++ )
		{
			for( int j = 0; j < iCols; j++ )
			{
				int iData = vvData[i][j];
				fwrite( &iData, sizeof(int), 1, infile );
			}
		}
	}
	else if( pcFormat[0] == 'f' )
	{
		for( int i = 0; i < iRows; i++ )
		{
			for( int j = 0; j < iCols; j++ )
			{
				double dData = vvData[i][j];
				fwrite( &dData, sizeof(double), 1, infile );
			}
		}
	}
	fclose( infile );

	return 0;
}

int
main_runMatching_tmi2017(
				 char *pcFileQueryImages,
				 char *pcFileBogusImages,
				 char *pcDirSaveResult,
				 char **pcFilesTrainingImages,
				 int iTrainingImageCount
				 )
{
	printf( "in\n" );

    string matcherType = defaultMatcherType;
    string fileWithQueryImages = defaultQueryImageName;
    string fileWithTrainImages = defaultFileWithTrainImages;
    string fileWithBogusImages = "";
    string dirToSaveResImages = defaultDirToSaveResImages;

	//matcherType = argv[3];
    fileWithQueryImages = pcFileQueryImages;
	fileWithBogusImages = pcFileBogusImages;
    dirToSaveResImages = pcDirSaveResult;
	fileWithTrainImages = pcFilesTrainingImages[0];

    Ptr<DescriptorMatcher> descriptorMatcher = DescriptorMatcher::create( matcherType );
    Ptr<DescriptorMatcher> descriptorMatcherBogus = DescriptorMatcher::create( matcherType );

 //   Mat maskImage = imread( fileMaskImage, CV_LOAD_IMAGE_GRAYSCALE);
	//if( maskImage.empty() )
	//{
	//	cout << "Mask image can not be read:" << fileMaskImage << endl << ">" << endl;
 //       return false;
	//}

	printf( "REading file:  %s\n", &fileWithQueryImages[0] );

	Mat queryImage;
	vector<Mat> queryImages;
    vector<string> queryImagesNames;
	vector< vector<float> > queryImagesTransform;
	vector<vector<KeyPoint> > queryKeypoints;
	vector<Mat> queryDescriptors;
	_readDescriptorDatabase( &fileWithQueryImages[0], queryImages, queryImagesNames,
						queryImagesTransform, queryKeypoints, queryDescriptors );

    vector<Mat> trainImages;
    vector<string> trainImagesNames;
	vector< vector<float> > trainImagesTransform;
	vector<vector<KeyPoint> > trainKeypoints;
	vector<Mat> trainDescriptors;
	_readDescriptorDatabase( &fileWithTrainImages[0], trainImages, trainImagesNames, trainImagesTransform, trainKeypoints, trainDescriptors );
	//_readDescriptorDatabase_range( &fileWithTrainImages[0], trainImages, trainImagesNames, trainImagesTransform, trainKeypoints, trainDescriptors, 320, 500);


	descriptorMatcher->add( trainDescriptors );


	vector<Mat> bogusImages;
    vector<string> bogusImagesNames;
	vector< vector<float> > bogusImagesTransform;
	vector<vector<KeyPoint> > bogusKeypoints;
	vector<Mat> bogusDescriptors;
	_readDescriptorDatabase( &fileWithBogusImages[0], bogusImages, bogusImagesNames,
						bogusImagesTransform, bogusKeypoints, bogusDescriptors );
	descriptorMatcherBogus->add( bogusDescriptors );

#ifdef WIN32
		_timeb t1, t2;
		_ftime( &t1 );
#endif

#ifdef WIN32
		_ftime( &t2 );
		 int t_diff = (int) (1000.0 * (t2.time - t1.time) + (t2.millitm - t1.millitm));   
		printf( "Training time: %d\n", t_diff );
#endif
	vector< vector<float> > vecOutputX;
	vector< vector<float> > vecOutputY;

	// Open & flush output
	string resultName1 = dirToSaveResImages + "\\matches1.txt";
	string resultName2 = dirToSaveResImages + "\\matches2.txt";
	FILE *outfile1;
	FILE *outfile2;
	outfile1 = fopen( resultName1.c_str(), "wt" );
	outfile2 = fopen( resultName2.c_str(), "wt" );
	if( !outfile1 )
		return -1;
	if( !outfile2 )
		return -1;
	fclose( outfile1 );
	fclose( outfile2 );

	// Read each image, match to training
	//vector<KeyPoint> queryKeypoints;
	//Mat queryDescriptors;
    vector<DMatch> matches;
    vector<vector<DMatch> > mmatches;
    vector<vector<DMatch> > mmatchesBogus;
	vector<int> vecImgMatches;
	vecImgMatches.resize( trainImages.size() );
	FILE *outfileTransform = fopen( "match_transform.txt", "wt" );

	//
	// *** November 2014, visualize correspondences here
	// 
	Mat mtVisualizeCorr( queryImages.size(), trainImages.size(), CV_32FC1 );
	Mat mtVisualizeCorrSelected( queryImages.size(), trainImages.size(), CV_32FC1 );
	Mat mtVisualizeCorrNorm( queryImages.size(), trainImages.size(), CV_32FC1 );
	memset( mtVisualizeCorr.data, 0, sizeof(float)*queryImages.size()*trainImages.size() );


	for( int i = 0; i < queryImages.size(); i++ )
	{
		for( int j = 0; j < vecImgMatches.size(); j++ )
			vecImgMatches[j] = 0;

		//queryKeypoints.clear();
		mmatches.clear();
		mmatchesBogus.clear();

#ifdef WIN32
		_ftime( &t1 );
#endif
		vector<KeyPoint> &queryKeypointsImg = queryKeypoints[i];
		Mat &queryDescriptorsImg = queryDescriptors[i];
		
		if( queryKeypointsImg.size() == 0 )
		{
			continue;
		}

		//featureDetector->detect( queryImages[i], queryKeypoints, maskImage );
		//descriptorExtractor->compute( queryImages[i], queryKeypoints, queryDescriptors );
		matchDescriptorsKNN( queryDescriptorsImg, mmatches, descriptorMatcher, 3 );
		matchDescriptorsKNN( queryDescriptorsImg, mmatchesBogus, descriptorMatcherBogus, 1 );

		// Remove matches to self, same image
		for( int j = 0; j < mmatches.size(); j++ )
		{
			for( int k = 1; k < mmatches[j].size(); k++ )
			//for( int k = 1; k < mmatches[k].size(); k++ )
			{
				if( mmatches[j][0].imgIdx == i )
				{
					if(  mmatches[j][k].imgIdx != i )
					{
						mmatches[j][0] =  mmatches[j][k];
						break;
					}
				}
			}
			// Disable for now
			if( 0 && mmatches[j][0].imgIdx == i )
			{
				mmatches[j][0].queryIdx = -1;
				mmatches[j][0].trainIdx = -1;
				mmatches[j][0].imgIdx = -1;
			}
		}

		vector< DMatch > vmMatches;
		//float fRatioThreshold = 0.90;
		//float fRatioThreshold = 0.95;
		float fRatioThreshold = 1.0;
		for( int j = 0; j < mmatches.size(); j++ )
		{
			float fRatio = mmatches[j][0].distance / mmatchesBogus[j][0].distance;
			//float fRatio = mmatches[j][0].distance / mmatches[j][1].distance;
			if( fRatio > fRatioThreshold )
			{
				mmatches[j][0].queryIdx = -1;
				mmatches[j][0].trainIdx = -1;
				mmatches[j][0].imgIdx = -1;
			}
			else
			{
				vmMatches.push_back( mmatches[j][0] );
			}
		}

#ifdef WIN32
		_ftime( &t2 );
		 int t_diff = (int) (1000.0 * (t2.time - t1.time) + (t2.millitm - t1.millitm));   
		printf( "Testing Time: %d\n", t_diff );
#endif
		// Do hough transform to prune
		//int iMatchCount = houghTransform( queryKeypointsImg, trainKeypoints, vmMatches, 615, 400 );
		//int iMatchCount = houghTransform( queryKeypointsImg, trainKeypoints, vmMatches, 385, 153 ); // AMIGO December 2012 - MICCAI 2013 results
		//int iMatchCount = houghTransform( queryKeypointsImg, trainKeypoints, vmMatches, 410, 300 ); // AMIGO January 2013
		//int iMatchCount = houghTransform( queryKeypointsImg, trainKeypoints, vmMatches, 200, 400 ); // Tamas phantom 2013

		int iMatchCount = houghTransform( queryKeypointsImg, trainKeypoints, vmMatches, 315, 300 ); // MNI

		// Tally votes, find frame with the most matches
		for( int j = 0; j < vmMatches.size(); j++ )
		{
			int iImg =  vmMatches[j].imgIdx;
			if( iImg >= 0 )
			{
				vecImgMatches[iImg]++;
			}
		}

		// Save image of counts
		for( int j = 0; j < vecImgMatches.size(); j++ )
		{
			int iCount = vecImgMatches[j];
			((float*)mtVisualizeCorr.data)[i*trainImages.size()+j] = iCount;
		}

		// 
		// OK - what precisely is going on in this scenario
		// Start a new function 2014 abover to figure this out... 
		//
		vector< int > vecImgMatchesSmooth;
		vecImgMatchesSmooth.resize( vecImgMatches.size(), 0 );

		int iMaxIndex = -1;
		int iMaxCount = -1;
		for( int j = 0; j < vecImgMatches.size(); j++ )
		{
			int iCount = vecImgMatches[j];
			if( iCount > iMaxCount )
			{
				iMaxCount = iCount;
				iMaxIndex = j;
			}
			vecImgMatchesSmooth[j] = vecImgMatches[j];
			if( j > 0 ) vecImgMatchesSmooth[j] += vecImgMatches[j-1];
			if( j < vecImgMatches.size()-1 ) vecImgMatchesSmooth[j] += vecImgMatches[j+1];
		}

		for( int j = 0; j < vecImgMatchesSmooth.size(); j++ )
		{
			vecImgMatches[j] = 0;
		}
		for( int j = 0; j < vecImgMatchesSmooth.size(); j++ )
		{
			if( vecImgMatchesSmooth[j] >= 2 )
			{
				// flag neighborhood
				vecImgMatches[j] = 1;
				if( j > 0 ) vecImgMatches[j-1]=1;
				if( j < vecImgMatches.size()-1 ) vecImgMatches[j+1]=1;
			}
		}

		// Save all matches
		vector< DMatch > vmMatchesSmooth;
		vmMatchesSmooth.clear();
		for( int j = 0; j < vecImgMatches.size(); j++ )
		{
			if( vecImgMatches[j] > 0 )
			{
				for( int k = 0; k < vmMatches.size(); k++ )
				{
					int iImg =  vmMatches[k].imgIdx;
					if( iImg == j )
					{
						vmMatchesSmooth.push_back( vmMatches[k] );
					}
				}
			}
		}

		if( vmMatchesSmooth.size() < 2 )
		{
			continue;
		}


		// Save image of counts
		for( int j = 0; j < vecImgMatches.size(); j++ )
		{
			int iCount = vecImgMatches[j];
			((float*)mtVisualizeCorrSelected.data)[i*trainImages.size()+j] = iCount;
		}
		// saveResultImages( queryImages[i], queryKeypointsImg, trainImages, trainKeypoints, vmMatchesSmooth, trainImagesNames, dirToSaveResImages, i );

		// Open & append result
		outfile1 = fopen( resultName1.c_str(), "a+" );
		outfile2 = fopen( resultName2.c_str(), "a+" );
		if( !outfile1 || !outfile2 )
			return -1;
		saveResultData(
			outfile1, outfile2,
			queryImagesNames, queryKeypoints,
			trainImagesNames, trainKeypoints,
			vmMatchesSmooth, i );
		fclose( outfile1 );
		fclose( outfile2 );

		// Save resulting transform - not really used
		if( 0 )
		{
			outfileTransform = fopen( "match_transform.txt", "a+" );
			fprintf( outfileTransform, "A:\t" );
			for( int k = 0; k < queryImagesTransform[i].size(); k++ )
			{
				fprintf( outfileTransform, "%f\t", queryImagesTransform[i][k] );
			}
			fprintf( outfileTransform, "B:\t" );
			for( int k = 0; k < trainImagesTransform[iMaxIndex].size(); k++ )
			{
				fprintf( outfileTransform, "%f\t", trainImagesTransform[iMaxIndex][k] );
			}
			fprintf( outfileTransform, "\n" );
			fclose( outfileTransform );
		}

		//
		// Add to linear solver & test
		//addToOutput(
		//	vmMatchesSmooth,
		//	queryImagesTransform[i], queryKeypoints[i],
		//	trainImagesTransform, trainKeypoints,
		//	 vecOutputX, vecOutputY );
		// For testing result of linear solver (doesn't work well)
		//test_output( vmMatchesSmooth, queryImagesTransform[i], queryKeypoints[i], trainImagesTransform, trainKeypoints );

		addToOutputNonLinear(
			vmMatchesSmooth,
			queryImagesTransform[i], queryKeypoints[i],
			trainImagesTransform, trainKeypoints,
			 vecOutputX, vecOutputY );
		
	}

	char pcFileName[400];
	cv::normalize( mtVisualizeCorr, mtVisualizeCorrNorm, 255, 0, CV_MINMAX  );
	imwrite( "correspondence_visualization.png", mtVisualizeCorrNorm );
	sprintf( pcFileName, "%s.correspondence_visualization.png", pcFileQueryImages ); 
	imwrite( pcFileName, mtVisualizeCorrNorm );

	cv::normalize( mtVisualizeCorrSelected, mtVisualizeCorrNorm, 255, 0, CV_MINMAX  );
	imwrite( "correspondence_visualization_selected.png", mtVisualizeCorrNorm );
	sprintf( pcFileName, "%s.correspondence_visualization_selected.png", pcFileQueryImages ); 
	imwrite( pcFileName, mtVisualizeCorrNorm );

	return 0;

	//writeArray( "non-linear.matches.X", vecOutputX, "f" );
	//writeArray( "non-linear.matches.Y", vecOutputY, "f" );


	sprintf( pcFileName, "%s.non-linear.matches.X", pcFileQueryImages ); 
	writeArray( pcFileName, vecOutputX, "f" );

	sprintf( pcFileName, "%s.non-linear.matches.Y", pcFileQueryImages ); 
	writeArray( pcFileName, vecOutputY, "f" );

    return 0;
}



//
// main()
//
// Main function for feature extraction, feature matching
// TMI 2017
//
int main(int argc, char** argv)
{
	if( argc <= 1 )
	{
		printf( "To prepare data (feature extraction), arguments: <keypoint> <descriptor> <listing file> <mask>\n" );
		printf( "To match data, arguments: <new prep-file> <bogus prep-file> <outdir> <pre-op1 prep-file> <pre-op2 prep-file> ... \n" );
		printf( "To convert MHA to PNG files, arguments: <mha file> <output code> \n" );
		printf( "To convert PNG/list text file to MHA, arguments: <text file in> <mha file out> <image data directory> \n" );
		return 0;
	}


	// readTrainTransforms_mha( filename, dirName, trainFilenames, trainImages, trainImageData );

	//convertToPng_mha( argv[1], argv[2] );

	//
	// IN USE: read in one mha / text file sequence, prepare descriptors.
	//
#ifdef FEATURE_EXTRACT
	main_prepareDescriptorDatabase( argv[1], argv[2], argv[3], argv[4] );
#endif
    
	//
	// IN USE: Run matching code between two images.
	//
#ifdef FEATURE_MATCH
	main_runMatching_tmi2017( argv[1], argv[2], argv[3], argv + 4, 1 );
#endif

    return 1;

}

